// Original implementation of the geometric definition in doi:10.1002/mp.70144.
// No vendor references. Include this file in a .NET Framework ESAPI project.
using System;
using System.Collections.Generic;
using System.Linq;

namespace PamBam
{
    public sealed class Rect
    {
        public double X1 { get; private set; } public double Y1 { get; private set; }
        public double X2 { get; private set; } public double Y2 { get; private set; }
        public Rect(double x1, double y1, double x2, double y2)
        { X1=x1; Y1=y1; X2=x2; Y2=y2; }
    }
    public sealed class Layer
    {
        public double[] BoundariesMm { get; private set; }
        public double[] BankAMm { get; private set; }
        public double[] BankBMm { get; private set; }
        public string Axis { get; private set; }
        public Layer(double[] boundariesMm, double[] bankAMm, double[] bankBMm, string axis="X")
        {
            if(boundariesMm==null || bankAMm==null || bankBMm==null) throw new ArgumentException("Missing layer data");
            BoundariesMm=(double[])boundariesMm.Clone(); BankAMm=(double[])bankAMm.Clone();
            BankBMm=(double[])bankBMm.Clone(); Axis=axis;
        }
    }
    public sealed class BeamResult
    {
        public double Bam { get; internal set; } public double TotalMu { get; internal set; }
        public string TargetKey { get; internal set; } public string Sampling { get; internal set; }
        public double[] Am { get; internal set; } public double[] WeightsMu { get; internal set; }
        internal BeamResult() { }
    }
    public static class Metrics
    {
        private static double Finite(double n)
        { if(double.IsNaN(n)||double.IsInfinity(n)) throw new ArgumentException("Finite numeric input required"); return n; }

        private static List<Rect> Rectangles(IEnumerable<Rect> input)
        {
            if(input==null) throw new ArgumentException("Missing geometry is not a closed aperture");
            var list=new List<Rect>();
            foreach(var r in input)
            {
                if(r==null) throw new ArgumentException("Missing rectangle");
                Finite(r.X1); Finite(r.X2); Finite(r.Y1); Finite(r.Y2);
                if(r.X1>r.X2 || r.Y1>r.Y2) throw new ArgumentException("Reversed rectangle bounds");
                if(r.X2>r.X1 && r.Y2>r.Y1) list.Add(r);
            }
            return list;
        }

        public static double UnionArea(IEnumerable<Rect> rectangles)
        {
            var r=Rectangles(rectangles);
            var xs=r.SelectMany(t=>new[]{t.X1,t.X2}).Distinct().OrderBy(x=>x).ToArray();
            double area=0;
            for(int i=1;i<xs.Length;i++)
            {
                double x1=xs[i-1],x2=xs[i],length=0,end=double.NegativeInfinity;
                foreach(var s in r.Where(t=>t.X1<=x1 && t.X2>=x2).OrderBy(t=>t.Y1).ThenBy(t=>t.Y2))
                { length+=Math.Max(0,s.Y2-Math.Max(s.Y1,end)); end=Math.Max(end,s.Y2); }
                area+=(x2-x1)*length;
            }
            return Finite(area);
        }

        public static List<Rect> Intersect(IEnumerable<Rect> first,IEnumerable<Rect> second)
        {
            var a=Rectangles(first); var b=Rectangles(second); var result=new List<Rect>();
            foreach(var r in a) foreach(var s in b)
            {
                double x1=Math.Max(r.X1,s.X1),y1=Math.Max(r.Y1,s.Y1);
                double x2=Math.Min(r.X2,s.X2),y2=Math.Min(r.Y2,s.Y2);
                if(x2>x1 && y2>y1) result.Add(new Rect(x1,y1,x2,y2));
            }
            return result;
        }

        private static List<Rect> Strips(Layer layer)
        {
            if(layer==null) throw new ArgumentException("Every physical MLC layer must be supplied");
            var e=layer.BoundariesMm;var a=layer.BankAMm;var b=layer.BankBMm;
            if(a.Length==0 || b.Length!=a.Length || e.Length!=a.Length+1 || (layer.Axis!="X"&&layer.Axis!="Y"))
                throw new ArgumentException("Invalid leaf geometry");
            foreach(double v in e.Concat(a).Concat(b)) Finite(v);
            var result=new List<Rect>();
            for(int i=0;i<a.Length;i++)
            {
                if(e[i+1]<=e[i] || a[i]>b[i]) throw new ArgumentException("Invalid leaf boundaries or banks");
                if(a[i]==b[i]) continue;
                result.Add(layer.Axis=="X" ? new Rect(a[i],e[i],b[i],e[i+1]) : new Rect(e[i],a[i],e[i+1],b[i]));
            }
            return result;
        }
        private static List<Rect> Limit(List<Rect> opening,Rect jaws,Rect fixedLimits)
        {
            if(jaws!=null) opening=Intersect(opening,new[]{jaws});
            if(fixedLimits!=null) opening=Intersect(opening,new[]{fixedLimits});
            return opening;
        }
        public static List<Rect> SingleLayerAperture(Layer layer,Rect jaws=null,Rect fixedLimits=null)
        { return Limit(Strips(layer),jaws,fixedLimits); }
        public static List<Rect> DualLayerAperture(Layer proximal,Layer distal,Rect jaws=null,Rect fixedLimits=null)
        { return Limit(Intersect(Strips(proximal),Strips(distal)),jaws,fixedLimits); }

        public static double ApertureModulation(IEnumerable<Rect> target,IEnumerable<Rect> aperture)
        {
            var t=Rectangles(target);var a=Rectangles(aperture);double total=UnionArea(t);
            if(total<=0) throw new ArgumentException("Target projection must have positive area");
            return Math.Max(0,Math.Min(1,1-UnionArea(Intersect(t,a))/total));
        }

        public static BeamResult CalculateBam(IEnumerable<IEnumerable<Rect>> targets,
            IEnumerable<IEnumerable<Rect>> apertures,IEnumerable<double> weightsMu,
            string targetKey,string sampling="weighted-samples")
        {
            if(targets==null||apertures==null||weightsMu==null) throw new ArgumentException("Missing samples");
            var t=targets.ToArray();var a=apertures.ToArray();var w=weightsMu.ToArray();
            if(w.Length==0||t.Length!=w.Length||a.Length!=w.Length||string.IsNullOrWhiteSpace(targetKey)||string.IsNullOrWhiteSpace(sampling))
                throw new ArgumentException("Complete matched samples and target key required");
            foreach(double v in w) if(Finite(v)<0) throw new ArgumentException("Negative sample MU");
            double total=Finite(w.Sum());
            if(total<=0) throw new ArgumentException("Beam MU must be positive");
            var am=t.Select((value,i)=>ApertureModulation(value,a[i])).ToArray();
            return new BeamResult {Bam=Math.Max(0,Math.Min(1,am.Select((value,i)=>value*(w[i]/total)).Sum())),TotalMu=total,
                TargetKey=targetKey,Am=am,WeightsMu=w,Sampling=sampling};
        }

        public static BeamResult BamFromControlPoints(IEnumerable<IEnumerable<Rect>> targets,
            IEnumerable<IEnumerable<Rect>> apertures,IEnumerable<double> cumulativeWeights,double totalMu,string targetKey)
        {
            if(cumulativeWeights==null) throw new ArgumentException("Missing cumulative weights");
            var c=cumulativeWeights.ToArray();Finite(totalMu);foreach(var v in c) Finite(v);
            if(c.Length<2||c[0]!=0||c[c.Length-1]<=0||totalMu<=0) throw new ArgumentException("Invalid cumulative weights or MU");
            var w=new double[c.Length];
            for(int i=1;i<c.Length;i++)
            {
                if(c[i]<c[i-1]) throw new ArgumentException("Decreasing cumulative weights");
                double mu=totalMu*((c[i]-c[i-1])/c[c.Length-1]);w[i-1]+=mu/2;w[i]+=mu/2;
            }
            return CalculateBam(targets,apertures,w,targetKey,"control-point-endpoint-trapezoid");
        }

        public static double CalculatePam(IEnumerable<BeamResult> beamResults)
        {
            if(beamResults==null) throw new ArgumentException("Missing beams");
            var beams=beamResults.ToArray();
            if(beams.Length==0) throw new ArgumentException("At least one complete treatment beam required");
            foreach(var b in beams)
                if(b==null||Finite(b.Bam)<0||b.Bam>1||Finite(b.TotalMu)<=0||string.IsNullOrWhiteSpace(b.TargetKey)||b.TargetKey!=beams[0].TargetKey)
                    throw new ArgumentException("Complete BAM, positive MU and common target required");
            double total=Finite(beams.Sum(b=>b.TotalMu));
            return Math.Max(0,Math.Min(1,beams.Sum(b=>b.Bam*(b.TotalMu/total))));
        }
    }
}
