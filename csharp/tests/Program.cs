using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Text.Json;
using PamBam;

class Program
{
    static int count;
    static void Near(double expected, double actual)
    { if (double.IsNaN(actual) || Math.Abs(expected-actual)>1e-10) throw new Exception("Numerical mismatch"); count++; }
    static void Reject(Action action)
    { try { action(); } catch (ArgumentException) { count++; return; } throw new Exception("Invalid input accepted"); }
    static double[] Numbers(JsonElement value) { return value.EnumerateArray().Select(x=>x.GetDouble()).ToArray(); }
    static Rect Rectangle(JsonElement value) { var p=Numbers(value); return new Rect(p[0],p[1],p[2],p[3]); }
    static int Main()
    {
        var t = new[]{new Rect(-10,-10,10,10)};
        var half = new[]{new Rect(-10,-10,0,10)};
        var closed = new Rect[0];
        var a = new Layer(new double[]{-10,0,10}, new double[]{-10,-10}, new double[]{10,10});
        var b = new Layer(new double[]{-15,-5,5,15}, new double[]{0,-5,0}, new double[]{5,5,5});
        Near(150, Metrics.UnionArea(Metrics.DualLayerAperture(a,b)));
        Near(.625, Metrics.ApertureModulation(t,Metrics.DualLayerAperture(a,b)));
        Near(150, Metrics.UnionArea(Metrics.DualLayerAperture(b,a)));
        Near(250, Metrics.UnionArea(Metrics.SingleLayerAperture(new Layer(new double[]{-10,0,10},new double[]{-20,-10},new double[]{20,10}),new Rect(-15,-5,15,5))));
        Near(200, Metrics.UnionArea(Metrics.SingleLayerAperture(new Layer(new double[]{-10,10},new double[]{-20},new double[]{20},"Y"),null,new Rect(-5,-10,5,10))));
        Near(0, Metrics.ApertureModulation(t,t)); Near(1, Metrics.ApertureModulation(t,closed));
        var overlapping = new[]{new Rect(0,0,2,2),new Rect(1,0,3,2)};
        Near(6, Metrics.UnionArea(overlapping));
        Near(1.0/3, Metrics.ApertureModulation(overlapping,new[]{new Rect(0,0,2,2)}));
        var targets = new[]{t,t,t}; var apertures = new[]{t,half,closed};
        var r = Metrics.BamFromControlPoints(targets,apertures,new double[]{0,.25,1},100,"target");
        Near(.625,r.Bam); Near(100,r.TotalMu); Near(12.5,r.WeightsMu[0]); Near(50,r.WeightsMu[1]); Near(37.5,r.WeightsMu[2]);
        Near(.625,Metrics.BamFromControlPoints(targets,apertures,new double[]{0,25,100},100,"target").Bam);
        Near(0,Metrics.BamFromControlPoints(targets,new[]{closed,t,t},new double[]{0,0,1},100,"target").Bam);
        var beamA = Metrics.CalculateBam(new[]{t},new[]{half},new double[]{100},"target");
        var beamB = Metrics.CalculateBam(new[]{t},new[]{t},new double[]{300},"target");
        Near(.125,Metrics.CalculatePam(new[]{beamA,beamB}));
        var closedBeam=Metrics.CalculateBam(Enumerable.Repeat(t,9),Enumerable.Repeat(closed,9),Enumerable.Repeat(1.0,9),"target");
        Near(1,Metrics.CalculatePam(new[]{closedBeam}));
        Reject(()=>Metrics.CalculatePam(new BeamResult[0]));
        Reject(()=>Metrics.CalculatePam(new[]{beamA,Metrics.CalculateBam(new[]{t},new[]{t},new double[]{100},"different")}));
        Reject(()=>Metrics.ApertureModulation(closed,t));
        Reject(()=>Metrics.DualLayerAperture(a,null));
        Reject(()=>Metrics.UnionArea(new[]{new Rect(0,0,double.PositiveInfinity,1)}));
        Reject(()=>Metrics.SingleLayerAperture(new Layer(new double[]{0,0},new double[]{0},new double[]{1})));
        foreach(double mu in new[]{0,-1,double.NaN,double.PositiveInfinity})
            Reject(()=>Metrics.BamFromControlPoints(new[]{t,t},new[]{t,t},new double[]{0,1},mu,"target"));
        foreach(var w in new[]{new double[]{0,-1},new double[]{0,double.NaN},new double[]{1,2},new double[]{0,0}})
            Reject(()=>Metrics.BamFromControlPoints(new[]{t,t},new[]{t,t},w,100,"target"));
        using(var fixtures=JsonDocument.Parse(File.ReadAllText(Path.Combine(AppContext.BaseDirectory,"geometry_cases.json"))))
        foreach(var c in fixtures.RootElement.EnumerateArray())
        {
            var ls=c.GetProperty("layers").EnumerateArray().Select(x=>new Layer(Numbers(x.GetProperty("edges")),Numbers(x.GetProperty("a")),Numbers(x.GetProperty("b")),x.TryGetProperty("axis",out var axis)?axis.GetString():"X")).ToArray();
            Rect jaw=c.TryGetProperty("jaws",out var j)?Rectangle(j):null;
            var opening=ls.Length==1?Metrics.SingleLayerAperture(ls[0],jaw):Metrics.DualLayerAperture(ls[0],ls[1],jaw);
            Near(c.GetProperty("area").GetDouble(),Metrics.UnionArea(opening));
            Near(c.GetProperty("am").GetDouble(),Metrics.ApertureModulation(c.GetProperty("target").EnumerateArray().Select(Rectangle),opening));
        }
        Console.WriteLine("PASS " + count + " C# analytical assertions");
        return 0;
    }
}
