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

namespace RSMC
{
  public class FilterProcessor
  {
    private FilterOptions options;
    public FilterProcessor(FilterOptions options)
    {
      this.options = options;
    }

    public bool Process()
    {
      Console.Out.WriteLine("filter process started at {0}", DateTime.Now);
      Stopwatch watch = new Stopwatch();
      watch.Start();

      var rproc = new Process
      {
        StartInfo = new ProcessStartInfo
        {
          FileName = options.GetRCommand(),
          Arguments = string.Format("--vanilla -f {0}", options.TargetRFile),
          UseShellExecute = false,
          RedirectStandardOutput = true,
          CreateNoWindow = true
        }
      };

      Console.Out.WriteLine("running command : " + rproc.StartInfo.FileName + " " + rproc.StartInfo.Arguments);
      try
      {
        if (!rproc.Start())
        {
          Console.Out.WriteLine("R command cannot be started, check your parameters and ensure that R is available.");
        }
      }
      catch (Exception ex)
      {
        Console.Out.WriteLine("R command cannot be started : ", ex.Message);
        return false;
      }

      try
      {
        string line;
        while ((line = rproc.StandardOutput.ReadLine()) != null)
        {
          Console.Out.WriteLine(line);
        }
      }
      catch (Exception ex)
      {
        Console.Out.WriteLine("R command error : ", ex.Message);
        return false;
      }

      watch.Stop();
      Console.Out.WriteLine("filter process ended at {0}, cost {1}", DateTime.Now, watch.Elapsed);

      return true;
    }
  }
}
