// PAM / BAM 1.1 -- read-only single-file ESAPI plug-in, selectable MLC profiles.
// Run this .cs file from Eclipse with one external photon plan open.
// API lengths and mesh positions are mm in DICOM coordinates. No dose required.
// The selected native ROI mesh is projected at EVERY CP onto the isocenter plane.
// The 2/1/0.5 mm setting is BEV sampling, not resampling of the native ROI mesh.
// AM = blocked target projection fraction. BAM uses half each incremental MU
// interval at both endpoints (trapezoidal CP approximation). PAM uses beam MU,
// not Beam.WeightFactor. No leaf-motion/dose/transmission model is evaluated.
// Scope: HFS, conventional IEC geometry, static couch yaw; Millennium/HD120
// and Halcyon SX dual-layer profiles. Auto or explicit profile selection.
// HD120: 14 x 5 mm + 32 x 2.5 mm + 14 x 5 mm, -110 to +110 mm at isocenter.
// Millennium120: 10x10 + 40x5 + 10x10 mm. Halcyon SX: 29/28 pairs, 10 mm
// widths with 5 mm stagger, both layers intersected and fixed 280 mm field.
// Dual-layer native ordering MUST be explicitly selected from local knowledge;
// it is not inferred from 57 leaves. Check the mapping in local commissioning.
// Known conflicting models and incomplete native arrays are rejected.
// Native source positions independently check the assumed beam frame at each CP.
// Unknown model aliases must be reviewed locally before adding an exact alias.
// Native Eclipse execution/commissioning is pending. See accompanying README.
// Sources: https://github.com/Kiragroh/PAM-BAM ; DOI:10.1002/mp.70144.
// Varian API reference 16.1; compiled against ESAPI 18.0.1.261 locally.
// Target .NET Framework 4.8 / x64, C# 5 syntax.
// No vendor libraries, patient data, external services or additional source files.
//
// MIT License -- Copyright (c) 2026 Maximilian Grohmann
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

using System;
using System.Collections.Generic;
using System.Linq;
#if !PAMBAM_CORE_TEST
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Globalization;
using System.Threading;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Threading;
using VMS.TPS.Common.Model.API;
using VMS.TPS.Common.Model.Types;
#endif

namespace PamBamSimple
{
    public sealed class InputError : Exception
    {
        public InputError(string message) : base(message) { }
    }
    public struct Vec
    {
        public double X, Y, Z;
        public Vec(double x, double y, double z) { X=x; Y=y; Z=z; }
        public static Vec operator +(Vec a, Vec b) { return new Vec(a.X+b.X,a.Y+b.Y,a.Z+b.Z); }
        public static Vec operator -(Vec a, Vec b) { return new Vec(a.X-b.X,a.Y-b.Y,a.Z-b.Z); }
        public static Vec operator *(Vec a, double b) { return new Vec(a.X*b,a.Y*b,a.Z*b); }
        public double Dot(Vec b) { return X*b.X+Y*b.Y+Z*b.Z; }
        public double Length { get { return Math.Sqrt(Dot(this)); } }
    }
    public sealed class Mesh
    {
        public Vec[] Vertices;
        public int[] Triangles;
        public void Validate()
        {
            Calculation.Check(Vertices!=null && Vertices.Length>=4 && Triangles!=null &&
                Triangles.Length>=12 && Triangles.Length%3==0, "Empty or incomplete ROI surface mesh.");
            Calculation.Check(Vertices.Length<=2000000 && Triangles.Length<=6000000,
                "ROI surface is too large for this simple viewer.");
            foreach(var v in Vertices) { Calculation.Finite(v.X); Calculation.Finite(v.Y); Calculation.Finite(v.Z); }
            foreach(int i in Triangles) Calculation.Check(i>=0 && i<Vertices.Length,"Invalid mesh triangle index.");
        }
    }
    public sealed class Sample
    {
        public double X0, Y0, Step;
        public int Width, Height;
        public bool[] Target;
    }
    public sealed class MlcLayer
    {
        public int[] Indices;
        public double[] Edges;
    }
    public sealed class MlcProfile
    {
        public string Name;
        public int Kind;
        public MlcLayer[] Layers;
        public double[] FixedLimits;
        public void Validate(int leafCount)
        {
            Calculation.Check(Layers!=null && (Layers.Length==1 || Layers.Length==2),"One or two complete MLC layers required.");
            var indices=new List<int>();
            foreach(var layer in Layers)
            {
                Calculation.Check(layer!=null && layer.Indices!=null && layer.Edges!=null &&
                    layer.Indices.Length>0 && layer.Edges.Length==layer.Indices.Length+1,"Incomplete MLC strip profile.");
                for(int i=0;i<layer.Edges.Length;i++)
                {
                    Calculation.Finite(layer.Edges[i]);
                    if(i>0) Calculation.Check(layer.Edges[i]>layer.Edges[i-1],"MLC strip boundaries must increase.");
                }
                indices.AddRange(layer.Indices);
            }
            Calculation.Check(indices.OrderBy(i=>i).SequenceEqual(Enumerable.Range(0,leafCount)),
                "MLC profile must map every native leaf pair exactly once. Check profile and native bank count.");
        }
    }
    public static class Calculation
    {
        public static void Check(bool condition,string message)
        { if(!condition) throw new InputError(message); }
        public static double Finite(double value)
        { Check(!double.IsNaN(value) && !double.IsInfinity(value),"Non-finite geometry or weight."); return value; }

        public static Vec[] Frame(double gantry,double couch,double collimator)
        {
            double g=Finite(gantry)*Math.PI/180, t=-Finite(couch)*Math.PI/180,
                c=Finite(collimator)*Math.PI/180;
            var basis=new[]{new Vec(1,0,0),new Vec(0,1,0),new Vec(0,0,1)};
            for(int i=0;i<3;i++)
            {
                Vec a=basis[i];
                // HFS @ Rz(-couch) @ Ry(gantry) @ Rz(collimator).
                var b=new Vec(a.X*Math.Cos(c)-a.Y*Math.Sin(c),a.X*Math.Sin(c)+a.Y*Math.Cos(c),a.Z);
                var d=new Vec(b.X*Math.Cos(g)+b.Z*Math.Sin(g),b.Y,-b.X*Math.Sin(g)+b.Z*Math.Cos(g));
                var e=new Vec(d.X*Math.Cos(t)-d.Y*Math.Sin(t),d.X*Math.Sin(t)+d.Y*Math.Cos(t),d.Z);
                basis[i]=new Vec(e.X,-e.Z,e.Y);
            }
            return basis;
        }

        public static double[] EndpointWeights(double[] cumulative)
        {
            Check(cumulative!=null && cumulative.Length>=2,"At least two control points are required.");
            foreach(double w in cumulative) Finite(w);
            Check(Math.Abs(cumulative[0])<1e-8 && cumulative[cumulative.Length-1]>0,
                "Cumulative MU must start at zero and end above zero.");
            double final=cumulative[cumulative.Length-1];
            var weights=new double[cumulative.Length];
            for(int i=1;i<cumulative.Length;i++)
            {
                double delta=cumulative[i]-cumulative[i-1];
                Check(delta>=0,"Cumulative control-point weights decrease.");
                weights[i-1]+=delta/(2*final); weights[i]+=delta/(2*final);
            }
            return weights;
        }

        public static double WeightedMean(double[] values,double[] weights)
        {
            Check(values!=null && weights!=null && values.Length>0 && values.Length==weights.Length,
                "Incomplete BAM/PAM results.");
            double sum=0, total=0;
            for(int i=0;i<values.Length;i++)
            {
                double v=Finite(values[i]),w=Finite(weights[i]);
                Check(v>=0 && v<=1 && w>=0,"Invalid BAM/PAM value or MU.");
                sum+=v*w; total+=w;
            }
            Check(total>0 && !double.IsInfinity(total),"Positive total MU required.");
            return Finite(sum/total);
        }

        public static double[] Hd120Boundaries()
        {
            var edges=new double[61]; edges[0]=-110;
            for(int i=0;i<60;i++) edges[i+1]=edges[i]+(i<14 || i>=46 ? 5.0 : 2.5);
            return edges;
        }

        public static double[] Millennium120Boundaries()
        {
            var edges=new double[61]; edges[0]=-200;
            for(int i=0;i<60;i++) edges[i+1]=edges[i]+(i<10 || i>=50 ? 10.0 : 5.0);
            return edges;
        }

        static string ModelKey(string model)
        { return new string((model??"").Where(char.IsLetterOrDigit).ToArray()).ToUpperInvariant(); }

        public static bool IsHd120Model(string model)
        {
            // Exact semantic aliases, never infer the MLC from "120" or leaf count.
            string key=ModelKey(model);
            return key=="HD120" || key=="VARIANHD120" || key=="HIGHDEFINITION120" ||
                key=="VARIANHIGHDEFINITION120";
        }

        public static int ModelKind(string model)
        {
            if(IsHd120Model(model)) return 2;
            string key=ModelKey(model);
            if(new[]{"MILLENNIUM120","MILLENIUM120","VARIANMILLENNIUM120","VARIANMILLENIUM120"}.Contains(key)) return 1;
            if(new[]{"SX1","SX2","HALCYON","HALCYONSX1","HALCYONSX2","VARIANHALCYON","VARIANSX1","VARIANSX2"}.Contains(key)) return 3;
            return 0;
        }

        public static MlcProfile ResolveProfile(string model,int selection,int dualOrder)
        {
            int recognized=ModelKind(model),kind=selection==0 ? recognized : selection;
            Check(kind>=1 && kind<=3,"Unknown MLC model: explicitly select a matching, locally verified profile.");
            Check(selection==0 || recognized==0 || recognized==kind,
                "Selected profile conflicts with the native MLC model. Choose Auto or the matching profile.");
            if(kind==1 || kind==2)
                return new MlcProfile { Kind=kind,Name=kind==1 ? "Millennium 120 (SD)" : "HD120",
                    Layers=new[]{new MlcLayer { Indices=Enumerable.Range(0,60).ToArray(),
                        Edges=kind==1 ? Millennium120Boundaries() : Hd120Boundaries() }} };
            Check(dualOrder>=1 && dualOrder<=3,
                "Dual-layer: select the locally verified native leaf order. 57 pairs alone do not identify the order.");
            int[] layer29,layer28;
            if(dualOrder==1)
            { layer29=Enumerable.Range(0,29).Select(i=>2*i).ToArray();layer28=Enumerable.Range(0,28).Select(i=>2*i+1).ToArray(); }
            else if(dualOrder==2)
            { layer29=Enumerable.Range(0,29).ToArray();layer28=Enumerable.Range(29,28).ToArray(); }
            else
            { layer28=Enumerable.Range(0,28).ToArray();layer29=Enumerable.Range(28,29).ToArray(); }
            return new MlcProfile { Kind=3,Name="Halcyon SX / "+(dualOrder==1 ? "interleaved" : dualOrder==2 ? "29+28" : "28+29"),FixedLimits=new[]{-140.0,-140,140,140},
                Layers=new[]{
                    new MlcLayer { Indices=layer29,Edges=Enumerable.Range(0,30).Select(i=>-145.0+10*i).ToArray() },
                    new MlcLayer { Indices=layer28,Edges=Enumerable.Range(0,29).Select(i=>-140.0+10*i).ToArray() }} };
        }

        public static Sample Project(Mesh mesh,Vec iso,double sad,Vec[] frame,double step,Action pulse)
        {
            Check(Finite(step)>0 && Finite(sad)>0,"Positive grid spacing and SAD are required.");
            var x=new double[mesh.Vertices.Length]; var y=new double[x.Length];
            double minX=double.PositiveInfinity,minY=minX,maxX=double.NegativeInfinity,maxY=maxX;
            for(int i=0;i<x.Length;i++)
            {
                var v=mesh.Vertices[i]-iso;
                double depth=sad-v.Dot(frame[2]);
                Check(depth>1,"ROI is at or beyond the source plane.");
                x[i]=Finite(v.Dot(frame[0])*sad/depth); y[i]=Finite(v.Dot(frame[1])*sad/depth);
                minX=Math.Min(minX,x[i]); maxX=Math.Max(maxX,x[i]);
                minY=Math.Min(minY,y[i]); maxY=Math.Max(maxY,y[i]);
            }
            double left=Math.Floor(minX/step), bottom=Math.Floor(minY/step);
            double w=Math.Ceiling(maxX/step)-left, h=Math.Ceiling(maxY/step)-bottom;
            Check(w>0 && h>0 && w*h<=600000,"Target projection too large or empty; adjust BEV grid.");
            var result=new Sample { X0=(left+.5)*step,Y0=(bottom+.5)*step,Step=step,
                Width=(int)w,Height=(int)h,Target=new bool[(int)(w*h)] };
            for(int n=0;n<mesh.Triangles.Length;n+=3)
            {
                if(n%6144==0) pulse();
                int a=mesh.Triangles[n], b=mesh.Triangles[n+1], c=mesh.Triangles[n+2];
                FillTriangle(result,x[a],y[a],x[b],y[b],x[c],y[c]);
            }
            Check(result.Target.Any(v=>v),"No target samples at this resolution; choose a finer grid.");
            return result;
        }

        static double Edge(double ax,double ay,double bx,double by,double x,double y)
        { return (bx-ax)*(y-ay)-(by-ay)*(x-ax); }

        static void FillTriangle(Sample grid,double ax,double ay,double bx,double by,double cx,double cy)
        {
            double area=Edge(ax,ay,bx,by,cx,cy);
            if(Math.Abs(area)<1e-10) return;
            double sign=area>0 ? 1 : -1;
            int x1=Math.Max(0,(int)Math.Ceiling((Math.Min(ax,Math.Min(bx,cx))-grid.X0)/grid.Step-1e-9));
            int x2=Math.Min(grid.Width-1,(int)Math.Floor((Math.Max(ax,Math.Max(bx,cx))-grid.X0)/grid.Step+1e-9));
            int y1=Math.Max(0,(int)Math.Ceiling((Math.Min(ay,Math.Min(by,cy))-grid.Y0)/grid.Step-1e-9));
            int y2=Math.Min(grid.Height-1,(int)Math.Floor((Math.Max(ay,Math.Max(by,cy))-grid.Y0)/grid.Step+1e-9));
            for(int j=y1;j<=y2;j++) for(int i=x1;i<=x2;i++)
            {
                int index=j*grid.Width+i;
                if(grid.Target[index]) continue;
                double x=grid.X0+i*grid.Step,y=grid.Y0+j*grid.Step;
                if(sign*Edge(ax,ay,bx,by,x,y)>=-1e-8 && sign*Edge(bx,by,cx,cy,x,y)>=-1e-8 &&
                    sign*Edge(cx,cy,ax,ay,x,y)>=-1e-8) grid.Target[index]=true;
            }
        }

        public static double BlockedFraction(Sample grid,float[,] leaves,double[] jaws)
        { return BlockedFraction(grid,leaves,jaws,ResolveProfile("HD120",2,0)); }

        public static double BlockedFraction(Sample grid,float[,] leaves,double[] jaws,MlcProfile profile)
        {
            Check(leaves!=null && leaves.GetLength(0)==2 && profile!=null,"Two native leaf banks and an explicit MLC profile are required.");
            profile.Validate(leaves.GetLength(1));
            jaws=profile.FixedLimits??jaws;
            Check(jaws!=null && jaws.Length==4,"Complete native jaws required.");
            foreach(double p in jaws) Finite(p);
            // ESAPI VRect order stored here: X1,Y1,X2,Y2.
            Check(jaws[0]<=jaws[2] && jaws[1]<=jaws[3],"Invalid signed jaw limits.");
            for(int i=0;i<leaves.GetLength(1);i++)
                Check(Finite(leaves[0,i])<=Finite(leaves[1,i]),"Crossed opposing leaf tips are unsupported.");
            int total=0,opened=0;
            for(int j=0;j<grid.Height;j++)
            {
                double y=grid.Y0+j*grid.Step;
                double left=jaws[0],right=jaws[2];
                bool yOpen=y>=jaws[1] && y<jaws[3];
                foreach(var layer in profile.Layers)
                {
                    int strip=Array.BinarySearch(layer.Edges,y);
                    if(strip<0) strip=~strip-1;
                    if(strip<0 || strip>=layer.Indices.Length) { yOpen=false;break; }
                    int leaf=layer.Indices[strip];
                    // Intersect all layer openings; never add or average areas.
                    left=Math.Max(left,leaves[0,leaf]);right=Math.Min(right,leaves[1,leaf]);
                }
                for(int i=0;i<grid.Width;i++)
                {
                    if(!grid.Target[j*grid.Width+i]) continue;
                    total++;
                    double x=grid.X0+i*grid.Step;
                    if(yOpen && x>=left && x<right) opened++;
                }
            }
            Check(total>0,"Target projection is empty.");
            return 1.0-(double)opened/total;
        }
    }

#if !PAMBAM_CORE_TEST
    public sealed class BeamRow
    {
        public string Beam { get; set; }
        public string MLC { get; set; }
        public string Profile { get; set; }
        public string MU { get; set; }
        public string CPs { get; set; }
        public string BAM { get; set; }
        public string Seconds { get; set; }
        public string Status { get; set; }
    }

    public sealed class Viewer
    {
        readonly ExternalPlanSetup plan;
        readonly Window window;
        readonly ComboBox roiBox=new ComboBox(),gridBox=new ComboBox(),mlcBox=new ComboBox(),orderBox=new ComboBox();
        readonly Button calculate=new Button { Content="Calculate" },cancel=new Button { Content="Cancel",IsEnabled=false };
        readonly TextBlock pam=new TextBlock { Text="Plan PAM: —",FontSize=22,FontWeight=FontWeights.Bold };
        readonly TextBlock status=new TextBlock { Text="Choose the target ROI and press Calculate.",TextWrapping=TextWrapping.Wrap };
        readonly TextBlock detail=new TextBlock { TextWrapping=TextWrapping.Wrap };
        readonly DataGrid table=new DataGrid { AutoGenerateColumns=false,IsReadOnly=true,CanUserAddRows=false };
        readonly ObservableCollection<BeamRow> rows=new ObservableCollection<BeamRow>();
        bool running, cancelled;

        public Viewer(ExternalPlanSetup currentPlan,Window host)
        {
            plan=currentPlan; window=host;
            window.Title="PAM / BAM 1.1 — ESAPI"; window.Width=1180; window.Height=690;
            window.MinWidth=880; window.MinHeight=510;
            var panel=new DockPanel { Margin=new Thickness(18) };
            var header=new StackPanel(); DockPanel.SetDock(header,Dock.Top); panel.Children.Add(header);
            header.Children.Add(new TextBlock { Text="Plan Aperture Modulation",FontSize=25,FontWeight=FontWeights.Bold });
            header.Children.Add(new TextBlock { Text="Read-only | SD / HD / Halcyon dual-layer | numerical approximation",Margin=new Thickness(0,4,0,15) });
            var controls=new StackPanel { Orientation=Orientation.Horizontal,Margin=new Thickness(0,0,0,12) };
            header.Children.Add(controls);
            controls.Children.Add(new TextBlock { Text="Target ROI",VerticalAlignment=VerticalAlignment.Center });
            roiBox.Width=260; roiBox.Margin=new Thickness(8,0,16,0);
            var targets=plan.StructureSet.Structures.Where(s=>!s.IsEmpty && s.HasSegment).OrderBy(s=>s.Id).ToArray();
            Calculation.Check(targets.Length>0,"No non-empty segmented ROIs in the planning structure set.");
            roiBox.ItemsSource=targets; roiBox.DisplayMemberPath="Id";
            roiBox.SelectedItem=targets.FirstOrDefault(s=>s.DicomType=="PTV")??targets[0];
            controls.Children.Add(roiBox);
            controls.Children.Add(new TextBlock { Text="BEV grid (mm)",VerticalAlignment=VerticalAlignment.Center });
            gridBox.ItemsSource=new[]{"2.0","1.0","0.5"}; gridBox.SelectedIndex=0;
            gridBox.Width=65;gridBox.Margin=new Thickness(8,0,14,0);controls.Children.Add(gridBox);
            calculate.Padding=new Thickness(12,3,12,3);cancel.Padding=new Thickness(12,3,12,3);
            cancel.Margin=new Thickness(8,0,0,0);controls.Children.Add(calculate);controls.Children.Add(cancel);
            var mlcControls=new StackPanel { Orientation=Orientation.Horizontal,Margin=new Thickness(0,0,0,12) };
            header.Children.Add(mlcControls);
            mlcControls.Children.Add(new TextBlock { Text="MLC profile",VerticalAlignment=VerticalAlignment.Center });
            mlcBox.ItemsSource=new[]{"Auto (native model)","TrueBeam SD / Millennium 120","TrueBeam HD / HD120","Dual-layer / Halcyon SX"};
            mlcBox.SelectedIndex=0;mlcBox.Width=255;mlcBox.Margin=new Thickness(8,0,16,0);mlcControls.Children.Add(mlcBox);
            mlcControls.Children.Add(new TextBlock { Text="Dual-layer leaf order",VerticalAlignment=VerticalAlignment.Center });
            orderBox.ItemsSource=new[]{"Select verified native order","Interleaved: 29 even / 28 odd","Layer blocks: 29 then 28","Layer blocks: 28 then 29"};
            orderBox.SelectedIndex=0;orderBox.Width=245;orderBox.Margin=new Thickness(8,0,0,0);mlcControls.Children.Add(orderBox);
            var footer=new StackPanel { Margin=new Thickness(0,12,0,0) };
            DockPanel.SetDock(footer,Dock.Bottom);panel.Children.Add(footer);
            footer.Children.Add(pam);footer.Children.Add(status);footer.Children.Add(detail);
            footer.Children.Add(new TextBlock { Text="0 = target projection fully open   •   1 = fully blocked\nMU-weighted geometry, not dose coverage. Dual-layer uses the shared opening of both layers.\nDual-layer order is an explicit local mapping; select it only after comparison with the TPS leaf display.",
                TextWrapping=TextWrapping.Wrap,Margin=new Thickness(0,12,0,0) });
            string[] names={"Beam","MLC","Profile","MU","CPs","BAM","Seconds","Status"};
            int[] widths={90,140,150,75,50,75,65,350};
            for(int i=0;i<names.Length;i++) table.Columns.Add(new DataGridTextColumn {
                Header=names[i],Binding=new Binding(names[i]),Width=widths[i] });
            table.ItemsSource=rows; panel.Children.Add(table);window.Content=panel;
            calculate.Click+=(s,e)=>Calculate();cancel.Click+=(s,e)=>cancelled=true;
            table.SelectionChanged+=(s,e)=> { var row=table.SelectedItem as BeamRow; if(row!=null) detail.Text=row.Beam+" / "+row.Profile+": "+row.Status; };
            roiBox.SelectionChanged+=(s,e)=>Invalidate();gridBox.SelectionChanged+=(s,e)=>Invalidate();
            mlcBox.SelectionChanged+=(s,e)=> { UpdateOrderControl();Invalidate(); };
            orderBox.SelectionChanged+=(s,e)=>Invalidate();UpdateOrderControl();
            window.Closing+=(s,e)=> { if(running) { cancelled=true;e.Cancel=true; } };
        }

        void Invalidate()
        { if(!running) { rows.Clear();pam.Text="Plan PAM: —";status.Text="Selection changed. Press Calculate.";detail.Text=""; } }
        void UpdateOrderControl()
        { orderBox.IsEnabled=!running && (mlcBox.SelectedIndex==3 || (mlcBox.SelectedIndex==0 &&
            plan.Beams.Any(b=>!b.IsSetupField && b.MLC!=null && Calculation.ModelKind(b.MLC.Model)==3))); }
        void Pulse()
        {
            // All API access stays on the owning STA; no worker, task or async API use.
            var frame=new DispatcherFrame();
            window.Dispatcher.BeginInvoke(DispatcherPriority.Background,new Action(()=>frame.Continue=false));
            Dispatcher.PushFrame(frame);
            if(cancelled) throw new OperationCanceledException();
        }
        static Vec From(VVector v) { return new Vec(v.x,v.y,v.z); }

        static Mesh ReadMesh(Structure target)
        {
            var native=target.MeshGeometry;
            Calculation.Check(native!=null,"Selected target has no surface mesh.");
            var mesh=new Mesh { Vertices=native.Positions.Select(p=>new Vec(p.X,p.Y,p.Z)).ToArray(),
                Triangles=native.TriangleIndices.ToArray() };
            mesh.Validate(); return mesh;
        }

        double CalculateBeam(Beam beam,Mesh mesh,double spacing,MlcProfile profile,out int ncp)
        {
            Calculation.Check(beam.MLC!=null,"A treatment MLC is required.");
            Calculation.Check(beam.Applicator==null && !beam.Blocks.Any() && beam.Compensator==null && !beam.Wedges.Any(),
                "Applicators, custom blocks, compensators and wedges are unsupported.");
            Calculation.Check(string.IsNullOrEmpty(beam.MotionCompensationTechnique),"Motion compensation is unsupported.");
            Calculation.Check(beam.EnergyModeDisplayName.IndexOf("X",StringComparison.OrdinalIgnoreCase)>=0,
                "Only photon treatment beams are supported.");
            var cps=beam.ControlPoints.ToArray();ncp=cps.Length;
            var weights=Calculation.EndpointWeights(cps.Select(c=>c.MetersetWeight).ToArray());
            var am=new double[ncp];
            double sad=Calculation.Finite(beam.TreatmentUnit.SourceAxisDistance);
            Vec iso=From(beam.IsocenterPosition);
            double couch=Calculation.Finite(cps[0].PatientSupportAngle);
            // Verify the whole gantry plane, also for a static field. A single
            // source direction cannot detect rotation around that direction.
            foreach(double angle in new[]{0.0,90.0})
            {
                Vec expected=iso+Calculation.Frame(angle,couch,0)[2]*sad;
                Calculation.Check((From(beam.GetSourceLocation(angle))-expected).Length<.1,
                    "Native source geometry disagrees with the HFS IEC gantry plane.");
            }
            for(int i=0;i<ncp;i++)
            {
                var cp=cps[i];
                Calculation.Check(Math.Abs(cp.PatientSupportAngle-couch)<1e-6,"Couch motion during a beam is unsupported.");
                // NaN translations mean not recorded, as permitted by ESAPI. Reject changes/mixed validity.
                double[] first={cps[0].TableTopLateralPosition,cps[0].TableTopLongitudinalPosition,cps[0].TableTopVerticalPosition};
                double[] now={cp.TableTopLateralPosition,cp.TableTopLongitudinalPosition,cp.TableTopVerticalPosition};
                for(int k=0;k<3;k++) Calculation.Check((double.IsNaN(first[k]) && double.IsNaN(now[k])) ||
                    (!double.IsNaN(first[k]) && !double.IsNaN(now[k]) && Math.Abs(first[k]-now[k])<1e-5),
                    "Changing or inconsistent table positions are unsupported.");
                var frame=Calculation.Frame(cp.GantryAngle,couch,cp.CollimatorAngle);
                var nativeSource=From(beam.GetSourceLocation(cp.GantryAngle));
                Calculation.Check((nativeSource-(iso+frame[2]*sad)).Length<.1,
                    "Native source position disagrees with the HFS IEC frame (pitch/roll or unsupported geometry).");
                status.Text=string.Format(CultureInfo.InvariantCulture,"Beam {0} | CP {1}/{2} | {3:0.0} mm BEV grid",beam.Id,i+1,ncp,spacing);
                Pulse();
                if(weights[i]==0) continue;
                var target=Calculation.Project(mesh,iso,sad,frame,spacing,Pulse);
                double[] nativeJaws=null;
                if(profile.FixedLimits==null)
                { var jaws=cp.JawPositions;nativeJaws=new[]{jaws.X1,jaws.Y1,jaws.X2,jaws.Y2}; }
                am[i]=Calculation.BlockedFraction(target,cp.LeafPositions,nativeJaws,profile);
            }
            return Calculation.WeightedMean(am,weights);
        }

        void Calculate()
        {
            if(running) return;
            running=true;cancelled=false;calculate.IsEnabled=false;cancel.IsEnabled=true;
            roiBox.IsEnabled=gridBox.IsEnabled=mlcBox.IsEnabled=orderBox.IsEnabled=false;rows.Clear();detail.Text="";pam.Text="Plan PAM: calculating…";
            int failed=0,excluded=0,total=0;var bams=new List<double>();var mus=new List<double>();
            BeamRow active=null;
            var timer=Stopwatch.StartNew();
            try
            {
                Calculation.Check(Thread.CurrentThread.GetApartmentState()==ApartmentState.STA,"Use the owning ESAPI STA thread.");
                Calculation.Check(plan.TreatmentOrientation.ToString()=="HeadFirstSupine","Only Head First Supine plans are supported.");
                var target=(Structure)roiBox.SelectedItem;
                double spacing=double.Parse((string)gridBox.SelectedItem,CultureInfo.InvariantCulture);
                status.Text="Reading the native target surface…";Pulse();
                var mesh=ReadMesh(target);
                foreach(var beam in plan.Beams)
                {
                    active=new BeamRow { Beam=beam.Id,MLC=beam.MLC==null ? "—" : beam.MLC.Model,Profile="—",MU="—",CPs="—",BAM="—",Seconds="—",Status="Reading…" };
                    rows.Add(active);var clock=Stopwatch.StartNew();
                    try
                    {
                        Pulse();
                        if(beam.IsSetupField) { excluded++;active.Status="Setup field excluded";continue; }
                        Calculation.Check(beam.Meterset.Unit==DosimeterUnit.MU,"Beam meterset is not in MU.");
                        double mu=Calculation.Finite(beam.Meterset.Value);
                        Calculation.Check(mu>=0,"Negative beam MU.");active.MU=mu.ToString("0.00",CultureInfo.InvariantCulture);
                        if(mu==0) { excluded++;active.Status="Zero-MU field excluded";continue; }
                        total++;int ncp;
                        Calculation.Check(beam.MLC!=null,"A treatment MLC is required.");
                        var profile=Calculation.ResolveProfile(beam.MLC.Model,mlcBox.SelectedIndex,orderBox.SelectedIndex);
                        active.Profile=profile.Name;
                        double bam=CalculateBeam(beam,mesh,spacing,profile,out ncp);
                        active.CPs=ncp.ToString(CultureInfo.InvariantCulture);
                        active.BAM=bam.ToString("0.0000",CultureInfo.InvariantCulture);
                        active.Status="OK (sampled geometry)";bams.Add(bam);mus.Add(mu);
                    }
                    catch(OperationCanceledException) { active.Status="Cancelled";throw; }
                    catch(Exception ex)
                    { failed++;active.Status=SafeMessage(ex); }
                    finally
                    { active.Seconds=clock.Elapsed.TotalSeconds.ToString("0.0",CultureInfo.InvariantCulture);table.Items.Refresh();table.ScrollIntoView(active); }
                }
                Calculation.Check(failed==0 && total>0 && total==bams.Count,
                    string.Format("Plan PAM unavailable: {0} failed field(s), {1} complete field(s).",failed,bams.Count));
                double value=Calculation.WeightedMean(bams.ToArray(),mus.ToArray());
                pam.Text=string.Format(CultureInfo.InvariantCulture,"Plan PAM: {0:0.0000}  ({1:0.00}%)",value,100*value);
                status.Text=string.Format(CultureInfo.InvariantCulture,"{0} | {1} fields | {2:0.00} MU | {3:0.0} mm BEV | {4:0.0} s | {5} excluded",target.Id,total,mus.Sum(),spacing,timer.Elapsed.TotalSeconds,excluded);
            }
            catch(OperationCanceledException)
            { pam.Text="Plan PAM: unavailable";status.Text="Cancelled. Completed field rows remain visible; no partial plan PAM."; }
            catch(Exception ex)
            { pam.Text="Plan PAM: unavailable";status.Text=SafeMessage(ex); }
            finally
            { running=false;calculate.IsEnabled=true;cancel.IsEnabled=false;roiBox.IsEnabled=gridBox.IsEnabled=mlcBox.IsEnabled=true;UpdateOrderControl(); }
        }
        static string SafeMessage(Exception ex)
        { return ex is InputError ? ex.Message : "ESAPI read failed ("+ex.GetType().Name+"). Check the current plan/ROI."; }
    }
#endif
}

#if !PAMBAM_CORE_TEST
namespace VMS.TPS
{
    public class Script
    {
        public void Execute(ScriptContext context,Window window)
        {
            try
            {
                var plan=context.PlanSetup as ExternalPlanSetup;
                PamBamSimple.Calculation.Check(plan!=null && plan.StructureSet!=null,"Open one external treatment plan with a structure set in Eclipse.");
                new PamBamSimple.Viewer(plan,window);
            }
            catch(Exception ex)
            { MessageBox.Show(ex is PamBamSimple.InputError ? ex.Message : "Cannot read the Eclipse plan context.","PAM / BAM"); }
        }
    }
}
#endif
