using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CommandLine;
using CommandLine.Text;
using System.IO;
using System.Reflection;
using CQS.Genome.Pileup;
using RCPA.Seq;
using CQS.Commandline;

namespace RSMC
{
  public class FilterOptions : AbstractProgramOptions
  {
    public FilterOptions()
      : base()
    {
      this.SourceRFile = Path.ChangeExtension(Assembly.GetExecutingAssembly().Location, ".r");
    }

    public string SourceRFile { get; set; }

    [Option('v', "pvalue", MetaValue = "DOUBLE", DefaultValue = 0.01, HelpText = "pvalue used for significance test")]
    public double PValue { get; set; }

    [Option('q', "base_quality", MetaValue = "INT", DefaultValue = 20, HelpText = "Minimum base quality for mpileup result filter")]
    public int MinimumBaseQuality { get; set; }

    [Option('o', "output", MetaValue = "FILE", Required = true, HelpText = "Output file")]
    public string OutputFile { get; set; }

    [Option('d', "candidate_dir", MetaValue = "DIRECTORY", Required = true, HelpText = "Directory including candidate site files")]
    public string CandidatesDirectory { get; set; }

    public string TargetRFile
    {
      get
      {
        return new FileInfo(this.OutputFile + ".r").FullName.Replace("\\", "/");
      }
    }

    public string GetRCommand()
    {
      return this.Config.FindOrCreate("R", "R").Command;
    }

    public override bool PrepareOptions()
    {
      if (!File.Exists(this.SourceRFile))
      {
        ParsingErrors.Add(string.Format("file not exists : {0}", this.SourceRFile));
        return false;
      }

      if (!this.IsPileup)
      {
        if (!Directory.Exists(this.CandidatesDirectory))
        {
          ParsingErrors.Add(string.Format("directory not exists : {0}", this.CandidatesDirectory));
          return false;
        }

        if (Directory.GetFiles(this.CandidatesDirectory, "*.wsm").Length == 0)
        {
          ParsingErrors.Add(string.Format("there is no candidates in directory : {0}", this.CandidatesDirectory));
          return false;
        }
      }

      try
      {
        var lines = File.ReadAllLines(this.SourceRFile);
        using (StreamWriter sw = new StreamWriter(this.TargetRFile))
        {
          sw.WriteLine("setwd(\"{0}\")", Path.GetFullPath(this.CandidatesDirectory).Replace("\\", "/"));
          sw.WriteLine("minscore<-{0}", this.MinimumBaseQuality);
          sw.WriteLine("filename<-\"{0}\"", Path.GetFullPath(this.OutputFile).Replace("\\", "/"));
          sw.WriteLine("pvalue<-{0}", this.PValue);
          bool setwd = true, minscore = true, filename = true, pvalue = true;
          foreach (var line in lines)
          {
            if (setwd && line.StartsWith("setwd("))
            {
              setwd = false;
              continue;
            }
            if (minscore && line.StartsWith("minscore<-"))
            {
              minscore = false;
              continue;
            }
            if (filename && line.StartsWith("filename<-"))
            {
              filename = false;
              continue;
            }
            if (pvalue && line.StartsWith("pvalue<-"))
            {
              pvalue = false;
              continue;
            }
            sw.WriteLine(line);
          }
        }

        return true;
      }
      catch (Exception ex)
      {
        ParsingErrors.Add(string.Format("create R file {0} failed: {1}", this.TargetRFile, ex.Message));
        return false;
      }
    }
  }
}
