// Vendor-free synthetic tests; compile with /define:PAMBAM_CORE_TEST.
using System;
using System.Collections.Generic;
using System.Linq;
using PamBamSimple;

class CoreTests
{
    static int checks;
    static void Assert(bool ok,string message) { checks++;if(!ok) throw new Exception(message); }
    static void Near(double a,double b,string message) { Assert(Math.Abs(a-b)<1e-8,message+": "+a+" / "+b); }
    static void Reject(Action action,string message)
    { bool rejected=false;try { action(); } catch(InputError) { rejected=true; } Assert(rejected,message); }
    static Mesh Boxes(Vec[] lows,Vec[] highs)
    {
        var vertices=new List<Vec>();var triangles=new List<int>();
        for(int i=0;i<lows.Length;i++)
        {
            Vec a=lows[i],b=highs[i];int offset=vertices.Count;
            vertices.AddRange(new[]{new Vec(a.X,a.Y,a.Z),new Vec(b.X,a.Y,a.Z),new Vec(b.X,b.Y,a.Z),new Vec(a.X,b.Y,a.Z),
                new Vec(a.X,a.Y,b.Z),new Vec(b.X,a.Y,b.Z),new Vec(b.X,b.Y,b.Z),new Vec(a.X,b.Y,b.Z)});
            triangles.AddRange(new[]{0,1,2,0,2,3,4,5,6,4,6,7,0,1,5,0,5,4,1,2,6,1,6,5,2,3,7,2,7,6,3,0,4,3,4,7}.Select(n=>n+offset));
        }
        var mesh=new Mesh { Vertices=vertices.ToArray(),Triangles=triangles.ToArray() };mesh.Validate();return mesh;
    }
    static bool RayBox(Vec source,Vec direction,Vec lower,Vec upper)
    {
        double[] s={source.X,source.Y,source.Z},d={direction.X,direction.Y,direction.Z};
        double[] lo={lower.X,lower.Y,lower.Z},hi={upper.X,upper.Y,upper.Z};
        double near=0,far=double.PositiveInfinity;
        for(int k=0;k<3;k++)
        {
            if(Math.Abs(d[k])<1e-12) { if(s[k]<lo[k] || s[k]>hi[k]) return false; }
            else
            { double a=(lo[k]-s[k])/d[k],b=(hi[k]-s[k])/d[k];near=Math.Max(near,Math.Min(a,b));far=Math.Min(far,Math.Max(a,b)); }
        }
        return far>near;
    }
    static float[,] Leaves(float left,float right)
    { var result=new float[2,60];for(int i=0;i<60;i++) { result[0,i]=left;result[1,i]=right; }return result; }
    static void TestProjection()
    {
        var lows=new[]{new Vec(-18,-7,-13),new Vec(6,-4,8)};
        var highs=new[]{new Vec(-6,9,-1),new Vec(19,6,18)};
        var mesh=Boxes(lows,highs);
        double[][] angles={new[]{0.0,0,0},new[]{37.0,23,61},new[]{90.0,0,90},new[]{180.0,31,17},new[]{270.0,0,0}};
        foreach(var a in angles)
        {
            Vec[] frame=Calculation.Frame(a[0],a[1],a[2]);double sad=1000;
            var grid=Calculation.Project(mesh,new Vec(),sad,frame,1,()=>{});
            var source=frame[2]*sad;
            bool equal=true;
            for(int j=0;j<grid.Height;j++) for(int i=0;i<grid.Width;i++)
            {
                double x=grid.X0+i*grid.Step,y=grid.Y0+j*grid.Step;
                Vec d=frame[0]*x+frame[1]*y-source;
                bool expected=false;
                for(int b=0;b<lows.Length;b++) expected|=RayBox(source,d,lows[b],highs[b]);
                if(expected!=grid.Target[j*grid.Width+i]) equal=false;
            }
            Assert(equal,"Perspective triangle raster vs independent ray/box oracle at gantry "+a[0]);
        }
        // An explicit through-hole remains empty: four boxes form a square ring.
        var ring=Boxes(new[]{new Vec(-20,-1,-20),new Vec(10,-1,-20),new Vec(-10,-1,-20),new Vec(-10,-1,10)},
                       new[]{new Vec(-10,1,20),new Vec(20,1,20),new Vec(10,1,-10),new Vec(10,1,20)});
        var rg=Calculation.Project(ring,new Vec(),1000,Calculation.Frame(0,0,0),2,()=>{});
        bool hole=true;
        for(int j=0;j<rg.Height;j++) for(int i=0;i<rg.Width;i++)
            if(Math.Abs(rg.X0+i*2)<8 && Math.Abs(rg.Y0+j*2)<8 && rg.Target[j*rg.Width+i]) hole=false;
        Assert(hole,"Projected hole is not convex-filled");
        var cube=Boxes(new[]{new Vec(-10,-10,-10)},new[]{new Vec(10,10,10)});
        var cg=Calculation.Project(cube,new Vec(),100,Calculation.Frame(0,0,0),1,()=>{});
        Assert(cg.Target.Count(v=>v)==22*22,"Analytical perspective cube: 100/90 magnification");
        Reject(()=>Calculation.Project(cube,new Vec(),10,Calculation.Frame(0,0,0),1,()=>{}),"Source plane crossing rejected");
        bool cancelled=false;
        try { Calculation.Project(cube,new Vec(),1000,Calculation.Frame(0,0,0),1,()=>{throw new OperationCanceledException();}); }
        catch(OperationCanceledException) { cancelled=true; }
        Assert(cancelled,"Projection cancellation propagates");
    }
    static void TestApertures()
    {
        var cube=Boxes(new[]{new Vec(-10,-10,-10)},new[]{new Vec(10,10,10)});
        var grid=Calculation.Project(cube,new Vec(),1000,Calculation.Frame(0,0,0),1,()=>{});
        var jaws=new[]{-100.0,-100,100,100};
        Near(Calculation.BlockedFraction(grid,Leaves(-100,100),jaws),0,"Open aperture");
        Near(Calculation.BlockedFraction(grid,Leaves(0,0),jaws),1,"Closed aperture");
        Near(Calculation.BlockedFraction(grid,Leaves(0,100),jaws),.5,"Half MLC opening");
        Near(Calculation.BlockedFraction(grid,Leaves(-100,100),new[]{0.0,-100,100,100}),.5,"Asymmetric native jaw order");
        Reject(()=>Calculation.BlockedFraction(grid,new float[2,57],jaws),"Dual-layer array rejected");
        Reject(()=>Calculation.BlockedFraction(grid,Leaves(10,-10),jaws),"Crossed tips rejected");
        var line=new Sample { X0=.1,Y0=-109.75,Step=.5,Width=1,Height=440,Target=Enumerable.Repeat(true,440).ToArray() };
        var leaves=Leaves(0,0);leaves[0,0]=-1;leaves[1,0]=1;
        Near(1-Calculation.BlockedFraction(line,leaves,new[]{-2.0,-120,2,120}),10.0/440,"Outer HD120 leaf has 5 mm width");
        leaves=Leaves(0,0);leaves[0,14]=-1;leaves[1,14]=1;
        Near(1-Calculation.BlockedFraction(line,leaves,new[]{-2.0,-120,2,120}),5.0/440,"Central HD120 leaf has 2.5 mm width");
    }
    static void Main()
    {
        double[] edges=Calculation.Hd120Boundaries();
        Assert(edges.Length==61,"HD120 pair count");
        Near(edges[0],-110,"HD120 lower extent");Near(edges[14],-40,"HD120 lower transition");
        Near(edges[46],40,"HD120 upper transition");Near(edges[60],110,"HD120 upper extent");
        Assert(Calculation.IsHd120Model("HD120"),"HD120 exact model");
        Assert(!Calculation.IsHd120Model("Millennium120"),"Never infer from 120");
        Assert(!Calculation.IsHd120Model("Halcyon"),"No guessed Halcyon profile");
        double[] weights=Calculation.EndpointWeights(new[]{0.0,20,60,100});
        double[] expected={.1,.3,.4,.2};for(int i=0;i<4;i++) Near(weights[i],expected[i],"Endpoint weight "+i);
        Near(Calculation.WeightedMean(new[]{.10,.55,.30,.15},weights),.325,"Published repository synthetic example");
        Near(Calculation.WeightedMean(new[]{.1,.5},new[]{100.0,300}),.4,"Plan MU weighting");
        var zero=Calculation.EndpointWeights(new[]{0.0,0,1,1});
        Near(zero[0],0,"No MU before first interval");Near(zero[3],0,"No MU after last interval");
        Reject(()=>Calculation.EndpointWeights(new[]{.1,1.0}),"Nonzero first cumulative weight rejected");
        Reject(()=>Calculation.EndpointWeights(new[]{0.0,1,.5}),"Decreasing weights rejected");
        Reject(()=>Calculation.WeightedMean(new[]{0.0},new[]{0.0}),"Zero total MU rejected");
        Reject(()=>Calculation.WeightedMean(new[]{double.NaN},new[]{1.0}),"Nonfinite result rejected");
        Near((Calculation.Frame(0,0,0)[2]-new Vec(0,-1,0)).Length,0,"Gantry 0 source anterior");
        Near((Calculation.Frame(90,0,0)[2]-new Vec(1,0,0)).Length,0,"Gantry 90 source patient-left");
        Near((Calculation.Frame(90,90,0)[2]-new Vec(0,0,-1)).Length,0,"Inverse couch yaw");
        Near((Calculation.Frame(0,0,90)[0]-new Vec(0,0,1)).Length,0,"Collimator rotation");
        TestProjection();TestApertures();
        Console.WriteLine("PASS: "+checks+" synthetic ESAPI GUI core assertions (no vendor runtime).");
    }
}
