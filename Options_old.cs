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

namespace wsmdetector
{
  public enum ProgramType
  {
    bam,
    mpileup,
    console
  }

  public class Options
  {
    public Options()
    {
      this.IgnoreInsertionDeletion = true;
      this.IgnoreTerminalBase = false;
      this.IgnoreN = false;
      this.RFile = Path.ChangeExtension(Assembly.GetExecutingAssembly().Location, ".r");
    }

    public bool IgnoreInsertionDeletion { get; set; }

    public bool IgnoreTerminalBase { get; set; }

    public bool IgnoreN { get; set; }

    public string RFile { get; set; }

    [Option('t', "type", Required = true, HelpText = "Where to read/generate mpileup result file (bam/mpileup/console)")]
    public ProgramType From { get; set; }

    [OptionList('b', "bamFiles", Separator = ',', Required = false, HelpText = "Bam files used in samtools mpileup (bam1,bam2), only two bam files are allowed")]
    public IList<string> MpileupBamFiles { get; set; }

    [Option('g', "genomeFasta", MetaValue = "FILE", Required = false, HelpText = "Genome fasta file for samtools mpileup")]
    public string MpileupGenomeFile { get; set; }

    [Option('n', "mpileupMinimumReadQuality", MetaValue = "INT", DefaultValue = 20, HelpText = "Minimum mapQ of read for samtools mpileup")]
    public int MpileupMinimumReadQuality { get; set; }

    [Option('m', "mpileupFile", MetaValue = "FILE", Required = false, HelpText = "Samtools mpileup result file")]
    public string MpileupResultFile { get; set; }

    [Option('v', "pvalue", MetaValue = "DOUBLE", DefaultValue = 0.01, HelpText = "pvalue used for significance test")]
    public double PValue { get; set; }

    [Option('q', "minimumBaseQuality", MetaValue = "INT", DefaultValue = 20, HelpText = "Minimum base quality for mpileup result filter")]
    public int MinimumBaseQuality { get; set; }

    [Option('a', "minimumPercentage", DefaultValue = 0.1, HelpText = "Minimum percentage of minor allele at at least one sample")]
    public double MinimumPercentageOfMinorAllele { get; set; }

    [Option('d', "minimumReadDepth", MetaValue = "INT", DefaultValue = 10, HelpText = "Minimum read depth of base passed mapping quality filter in each sample")]
    public int MinimumReadDepth { get; set; }

    [Option('p', "filterPosition", DefaultValue = true, HelpText = "Filter mpileup result based on position bias by fisher exact test")]
    public bool FilterPosition { get; set; }

    [Option('s', "filterStrand", DefaultValue = true, HelpText = "Filter mpileup result based on strand bias by fisher exact test")]
    public bool FilterStrand { get; set; }

    [Option('c', "threadCount", MetaValue = "INT", DefaultValue = 1, HelpText = "Number of thread")]
    public int ThreadCount { get; set; }

    [OptionList('r', "chromosomeNames", Separator = ',', Required = false, HelpText = "Chromosome names (separted by ',')")]
    public IList<string> ChromosomeNames { get; set; }

    [Option('o', "outputDirectory", Required = true, HelpText = "Output directory")]
    public string OutputDirectory { get; set; }

    public string OutputFileDirectory
    {
      get
      {
        return this.OutputDirectory + "/candidates";
      }
    }

    public string TargetRFile { get; set; }

    public string TargetRResultFile { get; set; }

    public string TargetAnnovarInputFile { get; set; }

    [ParserState]
    public IParserState LastParserState { get; set; }

    [HelpOption]
    public string GetUsage()
    {
      return HelpText.AutoBuild(this, (HelpText current) => HelpText.DefaultParsingErrorsHandler(this, current));
    }

    public PileupItemParser GetPileupItemParser()
    {
      return new PileupItemParser(this.MinimumReadDepth, 0, this.IgnoreInsertionDeletion, this.IgnoreN, this.IgnoreTerminalBase);
    }

    public bool ParseOptions(string[] args, ExecuteConfig config)
    {
      try
      {
        if (!CommandLine.Parser.Default.ParseArguments(args, this))
        {
          return false;
        }
      }
      catch (Exception ex)
      {
        Console.Out.WriteLine(ex.Message);
        Console.Out.WriteLine(this.GetUsage());
        return false;
      }

      if (!CreateOutputDirectory())
      {
        Console.Out.WriteLine(this.GetUsage());
        return false;
      }


      if (!WriteRFile())
      {
        Console.Out.WriteLine(this.GetUsage());
        return false;
      }

      switch (this.From)
      {
        case ProgramType.mpileup:
          if (null == this.MpileupResultFile)
          {
            Console.Out.WriteLine("ERROR: Mpileup file not defined.");
            Console.Out.WriteLine(this.GetUsage());
            return false;
          }
          if (!File.Exists(this.MpileupResultFile))
          {
            Console.Out.WriteLine("ERROR: Mpileup file not exists {0}.", this.MpileupResultFile);
            Console.Out.WriteLine(this.GetUsage());
            return false;
          }
          Console.Out.WriteLine("#mpileup file: " + this.MpileupResultFile);
          break;
        case ProgramType.bam:
          if (!File.Exists(config.SamtoolsExecute))
          {
            Console.Out.WriteLine("ERROR: Samtools not exist: {0}", config.SamtoolsExecute);
            Console.Out.WriteLine(this.GetUsage());
            return false;
          }

          if (null == this.MpileupBamFiles || 0 == this.MpileupBamFiles.Count)
          {
            Console.Out.WriteLine("ERROR: Bam files not defined.");
            Console.Out.WriteLine(this.GetUsage());
            return false;
          }
          if (2 != this.MpileupBamFiles.Count)
          {
            Console.Out.WriteLine("ERROR: Must be two indexed bam files for samtools mpileup.");
            Console.Out.WriteLine(this.GetUsage());
            return false;
          }
          foreach (var bamfile in this.MpileupBamFiles)
          {
            if (!File.Exists(bamfile))
            {
              Console.Out.WriteLine("ERROR: Bam file is not exists {0}", bamfile);
              Console.Out.WriteLine(this.GetUsage());
              return false;
            }
          }
          if (null == this.MpileupGenomeFile)
          {
            Console.Out.WriteLine("ERROR: Genome fasta file not defined.");
            Console.Out.WriteLine(this.GetUsage());
            return false;
          }
          if (!File.Exists(this.MpileupGenomeFile))
          {
            Console.Out.WriteLine("ERROR: Genome fasta file not exists {0}.", this.MpileupGenomeFile);
            Console.Out.WriteLine(this.GetUsage());
            return false;
          }

          if (this.ThreadCount >= 2)
          {
            Console.WriteLine("Checking chromosome names for thread mode ...");
            if (this.ChromosomeNames == null || this.ChromosomeNames.Count == 0)
            {
              var fai = this.MpileupGenomeFile + ".fai";
              if (File.Exists(fai))
              {
                var lines = File.ReadAllLines(fai);
                this.ChromosomeNames = lines.ToList().ConvertAll(m =>
                {
                  var pos = m.IndexOfAny(new char[] { '\t', ' ' });
                  if (pos == -1)
                  {
                    return m;
                  }
                  else
                  {
                    return m.Substring(0, pos);
                  }
                });
              }
              else
              {
                Console.WriteLine("Reading chromosome names from fasta file ...");
                this.ChromosomeNames = SequenceUtils.ReadFastaNames(this.MpileupGenomeFile);
                if (this.ChromosomeNames.Count == 0)
                {
                  Console.WriteLine("ERROR: Genome fasta file doesn't contain chromosome names, {0}.", this.MpileupGenomeFile);
                  Console.Out.WriteLine(this.GetUsage());
                  return false;
                }
              }
            }

            foreach (var chr in this.ChromosomeNames)
            {
              Console.WriteLine(chr);
            }
          }
          else
          {
            if (this.ChromosomeNames != null && this.ChromosomeNames.Count > 0)
            {
              Console.Out.WriteLine("#mpileup chromosome names: " + this.ChromosomeNames.Merge(","));
            }
          }

          Console.Out.WriteLine("#mpileup bam files: " + this.MpileupBamFiles.Merge(","));
          Console.Out.WriteLine("#mpileup minimum read quality: " + this.MpileupMinimumReadQuality);
          Console.Out.WriteLine("#mpileup genome fasta: " + this.MpileupGenomeFile);

          break;
        case ProgramType.console:
          Console.Out.WriteLine("#mpileup from console.");
          break;
      }

      return true;
    }

    private bool CreateOutputDirectory()
    {
      if (!Directory.Exists(this.OutputDirectory))
      {
        try
        {
          Directory.CreateDirectory(this.OutputDirectory);
        }
        catch (Exception ex)
        {
          Console.Out.WriteLine("Cannot create directory {0} : {1}", this.OutputDirectory, ex.Message);
          return false;
        }
      }

      if (!Directory.Exists(this.OutputFileDirectory))
      {
        try
        {
          Directory.CreateDirectory(this.OutputFileDirectory);
        }
        catch (Exception ex)
        {
          Console.Out.WriteLine("Cannot create directory {0} : {1}", this.OutputFileDirectory, ex.Message);
          return false;
        }
      }

      return true;
    }

    private bool WriteRFile()
    {
      this.TargetRFile = new FileInfo(this.OutputDirectory + "/" + Path.GetFileName(this.RFile)).FullName.Replace("\\", "/");
      if (!File.Exists(this.RFile))
      {
        Console.Out.WriteLine("File not exists : {0}", this.RFile);
        return false;
      }

      try
      {
        this.TargetRResultFile = new FileInfo(Path.GetDirectoryName(this.OutputFileDirectory) + "\\" + Path.GetFileName(this.OutputDirectory) + ".tsv").FullName.Replace("\\", "/");
        this.TargetAnnovarInputFile = this.TargetRResultFile + ".annovia";

        var lines = File.ReadAllLines(this.RFile);
        using (StreamWriter sw = new StreamWriter(this.TargetRFile))
        {
          sw.WriteLine("setwd(\"" + Path.GetFullPath(this.OutputFileDirectory).Replace("\\", "/") + "\")");
          sw.WriteLine("minscore<-{0}", this.MinimumBaseQuality);
          sw.WriteLine("filename<-\"../{0}\"", Path.GetFileName(this.TargetRResultFile));
          sw.WriteLine("annovarfilename<-\"../{0}\"", Path.GetFileName(this.TargetAnnovarInputFile));
          sw.WriteLine("pvalue<-{0}", this.PValue);
          bool setwd = true, minscore = true, filename = true, pvalue = true, annovar = true;
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
            if (annovar && line.StartsWith("annovarfilename<-"))
            {
              annovar = false;
              continue;
            }
            sw.WriteLine(line);
          }
        }

        return true;
      }
      catch (Exception ex)
      {
        Console.Out.WriteLine("Create R file failed: ", ex.Message);
        return false;
      }
    }
  }
}
