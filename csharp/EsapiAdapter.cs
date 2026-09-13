// Compile only inside your ESAPI project; vendor assemblies are not distributed.
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using VMS.TPS.Common.Model.API;
using VMS.TPS.Common.Model.Types;

namespace PamBam
{
    public sealed class EsapiLayerProfile
    {
        public int[] SourceLeafIndices;
        public double[] BoundariesMm;
        public string Axis="X";
    }
    public sealed class EsapiProfile
    {
        public string ExactMlcModel;
        public EsapiLayerProfile[] Layers;
        // Explicitly None for jawless devices; Fixed limits are not physical jaws.
        public string JawMode;
    }
    public static class EsapiAdapter
    {
        public static BeamResult CalculateBeam(Beam beam, EsapiProfile profile,
            Func<Beam,ControlPoint,IEnumerable<Rect>> targetProjectionAtControlPoint, string targetKey)
        {
            if(Thread.CurrentThread.GetApartmentState()!=ApartmentState.STA)
                throw new InvalidOperationException("Use the owning ESAPI STA thread");
            if(beam==null||beam.IsSetupField||beam.MLC==null||profile==null||targetProjectionAtControlPoint==null)
                throw new ArgumentException("Treatment MLC beam, profile and target projector required");
            if(string.IsNullOrWhiteSpace(profile.ExactMlcModel)||beam.MLC.Model!=profile.ExactMlcModel)
                throw new ArgumentException("Native MLC model must match the explicit profile");
            if(beam.Blocks.Any()||beam.Applicator!=null)
                throw new ArgumentException("Blocking accessories are outside this adapter geometry contract");
            if(beam.Meterset.Unit!=DosimeterUnit.MU) throw new ArgumentException("MU meterset required");
            var cps=beam.ControlPoints.ToArray();
            var targets=new List<IEnumerable<Rect>>();var apertures=new List<IEnumerable<Rect>>();
            foreach(var cp in cps)
            {
                apertures.Add(Aperture(cp.LeafPositions,cp.JawPositions,profile));
                // Copy immediately; never retain vendor objects for background workers.
                var projection=targetProjectionAtControlPoint(beam,cp);
                if(projection==null) throw new ArgumentException("Missing target projection");
                targets.Add(projection.ToArray());
            }
            return Metrics.BamFromControlPoints(targets,apertures,cps.Select(c=>c.MetersetWeight),beam.Meterset.Value,targetKey);
        }

        public static List<Rect> Aperture(float[,] banks,VRect<double> nativeJaws,EsapiProfile profile)
        {
            if(banks==null||banks.GetLength(0)!=2||profile==null||profile.Layers==null||
                (profile.Layers.Length!=1&&profile.Layers.Length!=2)||
                profile.Layers.Any(l=>l==null||l.SourceLeafIndices==null||l.BoundariesMm==null))
                throw new ArgumentException("Complete native banks and explicit layer profiles required");
            var indices=profile.Layers.SelectMany(l=>l.SourceLeafIndices).OrderBy(i=>i).ToArray();
            if(!indices.SequenceEqual(Enumerable.Range(0,banks.GetLength(1))))
                throw new ArgumentException("Every native leaf pair must be assigned exactly once");
            var layers=profile.Layers.Select(l=>new Layer(l.BoundariesMm,
                l.SourceLeafIndices.Select(i=>(double)banks[0,i]).ToArray(),
                l.SourceLeafIndices.Select(i=>(double)banks[1,i]).ToArray(),l.Axis)).ToArray();
            if(profile.JawMode!="None"&&profile.JawMode!="Physical"&&profile.JawMode!="Fixed")
                throw new ArgumentException("Explicit jaw mode required");
            Rect jaws=null, limits=null;
            if(profile.JawMode=="Physical") jaws=new Rect(nativeJaws.X1,nativeJaws.Y1,nativeJaws.X2,nativeJaws.Y2);
            if(profile.JawMode=="Fixed") limits=new Rect(nativeJaws.X1,nativeJaws.Y1,nativeJaws.X2,nativeJaws.Y2);
            return layers.Length==1 ? Metrics.SingleLayerAperture(layers[0],jaws,limits)
                : Metrics.DualLayerAperture(layers[0],layers[1],jaws,limits);
        }
    }
}
