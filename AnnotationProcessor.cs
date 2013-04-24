using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections.Concurrent;
using CQS.Genome.Pileup;
using CQS.Genome.Statistics;
using System.Threading;
using System.Diagnostics;
using System.IO;
using CQS.Genome.Annotation;
using CQS;

namespace RSMC
{
  public class AnnotationProcessor
  {
    private AnnotationOptions options;
    public AnnotationProcessor(AnnotationOptions options)
    {
      this.options = options;
    }

    public bool Process()
    {
      Console.Out.WriteLine("annotation process started at {0}", DateTime.Now);
      Stopwatch watch = new Stopwatch();
      watch.Start();

      List<IAnnotationCsvExporter> exporters = new List<IAnnotationCsvExporter>();

      if (options.Annovar)
      {
        if (!Directory.Exists(options.AnnotationDirectory))
        {
          try
          {
            Directory.CreateDirectory(options.AnnotationDirectory);
          }
          catch (Exception ex)
          {
            Console.Out.WriteLine("Cannot create directory {0} : {1}", options.AnnotationDirectory, ex.Message);
            return false;
          }
        }

        using (StreamWriter sw = new StreamWriter(options.AnnovarInputFile))
        {
          using (StreamReader sr = new StreamReader(options.InputFile))
          {
            string line = sr.ReadLine();
            while ((line = sr.ReadLine()) != null)
            {
              var parts = line.Split(',');
              if (parts.Length > 4)
              {
                sw.WriteLine(parts.Take(5).Merge('\t'));
              }
            }
          }
        }

        var perlproc = new Process
        {
          StartInfo = new ProcessStartInfo
          {
            FileName = "summarize_annovar.pl",
            Arguments = options.AnnovarParameter,
            UseShellExecute = false,
            RedirectStandardOutput = true,
            CreateNoWindow = true
          }
        };

        Console.Out.WriteLine("running command : " + perlproc.StartInfo.FileName + " " + perlproc.StartInfo.Arguments);
        try
        {
          if (!perlproc.Start())
          {
            Console.Out.WriteLine("Annovar command cannot be started, check your parameters and ensure that annovar and perl are available.");
          }
        }
        catch (Exception ex)
        {
          Console.Out.WriteLine("Annovar command cannot be started : ", ex.Message);
          return false;
        }

        try
        {
          string line;
          while ((line = perlproc.StandardOutput.ReadLine()) != null)
          {
            Console.Out.WriteLine(line);
          }
        }
        catch (Exception ex)
        {
          Console.Out.WriteLine("Annovar command error : ", ex.Message);
          return false;
        }

        var annovarResultFile = options.AnnovarOutputFile + ".genome_summary.csv";
        if (!File.Exists(annovarResultFile))
        {
          Console.Out.WriteLine("Annovar might be failed: cannot find annovar result {0} ", annovarResultFile);
          return false;
        }

        exporters.Add(new AnnovarExporter(annovarResultFile, (m, n) => GetKey(m, n)));
      }

      if (options.Rnaediting)
      {
        exporters.Add(new RnaeditingExporter(options.RnaeditingDatabase, (m, n) => GetKey(m, n)));
      }

      if (options.Distance)
      {
        if (!string.IsNullOrEmpty(options.DistanceJunctionBed))
        {
          exporters.Add(new JunctionDistanceExporter(options.DistanceJunctionBed));
        }
        if (!string.IsNullOrEmpty(options.DistanceInsertionBed))
        {
          exporters.Add(new InsertionDistanceExporter(options.DistanceInsertionBed));
        }
        if (!string.IsNullOrEmpty(options.DistanceDeletionBed))
        {
          exporters.Add(new DeletionDistanceExporter(options.DistanceDeletionBed));
        }
      }

      var resultfile = Path.ChangeExtension(options.InputFile, ".annotation.csv");
      Console.WriteLine("writing merged result " + resultfile + " ...");
      using (StreamWriter sw = new StreamWriter(resultfile))
      {
        using (StreamReader sr = new StreamReader(options.InputFile))
        {
          string line = sr.ReadLine();
          var parts = line.Split('\t');

          sw.WriteLine("{0},{1}", parts.Merge(','), (from exporter in exporters
                                                     select exporter.GetHeader()).Merge(','));

          var delimiter = options.InputFile.EndsWith(".tsv") ? '\t' : ',';

          while ((line = sr.ReadLine()) != null)
          {
            parts = line.Split(delimiter);
            var chr = parts[0];
            var position = long.Parse(parts[1]);
            //Console.WriteLine("{0}\t{1}", chr, position);
            sw.WriteLine("{0},{1}", parts.Merge(','), (from exporter in exporters
                                                       select exporter.GetValue(chr, position, position)).Merge(','));
          }
        }
      }

      watch.Stop();
      Console.Out.WriteLine("annotation process ended at {0}, cost {1}", DateTime.Now, watch.Elapsed);

      return true;
    }

    private string GetKey(string chr, long position)
    {
      return string.Format("{0}_{1}", chr, position);
    }
  }
}
