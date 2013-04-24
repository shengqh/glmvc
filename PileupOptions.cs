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
using CQS;

namespace RSMC
{
  public class PileupOptions : AbstractProgramOptions
  {
    public enum DataSourceType
    {
      bam,
      mpileup,
      console
    }

    public PileupOptions()
    {
      this.IgnoreInsertionDeletion = true;
      this.IgnoreTerminalBase = false;
      this.IgnoreN = false;
    }

    public bool IgnoreInsertionDeletion { get; set; }

    public bool IgnoreTerminalBase { get; set; }

    public bool IgnoreN { get; set; }

    [Option('t', "type", MetaValue = "STRING", Required = true, HelpText = "Where to read/generate mpileup result file (bam/mpileup/console)")]
    public DataSourceType From { get; set; }

    [OptionList('b', "bam", MetaValue = "BAM1,BAM2", Separator = ',', Required = false, HelpText = "Bam files used in samtools mpileup (INPUT MUST BE normal_bam,tumor_bam), only two bam files are allowed and the second one will be assumed as from tumor sample which contains somatic mutation")]
    public IList<string> BamFiles { get; set; }

    [Option('f', "fasta", MetaValue = "FILE", Required = false, HelpText = "Genome fasta file for samtools mpileup")]
    public string GenomeFastaFile { get; set; }

    [Option('n', "read_quality", MetaValue = "INT", DefaultValue = 20, HelpText = "Minimum mapQ of read for samtools mpileup")]
    public int MpileupMinimumReadQuality { get; set; }

    [Option('m', "mpileup", MetaValue = "FILE", Required = false, HelpText = "Samtools mpileup result file")]
    public string MpileupFile { get; set; }

    [Option('e', "pvalue", MetaValue = "DOUBLE", DefaultValue = 0.01, HelpText = "pvalue used for significance test")]
    public double PValue { get; set; }

    [Option('q', "base_quality", MetaValue = "INT", DefaultValue = 20, HelpText = "Minimum base quality for mpileup result filter")]
    public int MinimumBaseQuality { get; set; }

    [Option('g', "percentage", MetaValue = "DOUBLE", DefaultValue = 0.1, HelpText = "Minimum percentage of minor allele at at tumor sample (second sample)")]
    public double MinimumPercentageOfMinorAllele { get; set; }

    [Option('d', "read_depth", MetaValue = "INT", DefaultValue = 10, HelpText = "Minimum read depth of base passed mapping quality filter in each sample")]
    public int MinimumReadDepth { get; set; }

    [Option('p', "filter_position", DefaultValue = true, HelpText = "Filter mpileup result based on position bias by fisher exact test")]
    public bool FilterPosition { get; set; }

    [Option('s', "filter_strand", DefaultValue = true, HelpText = "Filter mpileup result based on strand bias by fisher exact test")]
    public bool FilterStrand { get; set; }

    [Option('c', "thread_count", MetaValue = "INT", DefaultValue = 1, HelpText = "Number of thread")]
    public int ThreadCount { get; set; }

    [OptionList('r', "chromosomes", MetaValue = "STRING", Separator = ',', Required = false, HelpText = "Chromosome names (separted by ',')")]
    public IList<string> ChromosomeNames { get; set; }

    [Option('o', "output", MetaValue = "DIRECTORY", Required = true, HelpText = "Output directory")]
    public string OutputDirectory { get; set; }

    public string CandidatesDirectory
    {
      get
      {
        return this.OutputDirectory + "/candidates";
      }
    }

    public PileupItemParser GetPileupItemParser()
    {
      return new PileupItemParser(this.MinimumReadDepth, this.MinimumBaseQuality, this.IgnoreInsertionDeletion, this.IgnoreN, this.IgnoreTerminalBase);
    }

    public string GetSamtoolsCommand()
    {
      return this.Config.FindOrCreate("samtools", "samtools").Command;
    }

    public override bool PrepareOptions()
    {
      if (!PrepareOutputDirectory())
      {
        return false;
      }

      switch (this.From)
      {
        case DataSourceType.mpileup:
          if (null == this.MpileupFile)
          {
            ParsingErrors.Add("Mpileup file not defined.");
            return false;
          }
          if (!File.Exists(this.MpileupFile))
          {
            ParsingErrors.Add(string.Format("Mpileup file not exists {0}.", this.MpileupFile));
            return false;
          }
          Console.Out.WriteLine("#mpileup file: " + this.MpileupFile);
          break;
        case DataSourceType.bam:
          if (!Config.Programs.ContainsKey("samtools"))
          {
            Config.Programs["samtools"] = new ProgramConfig("samtools", "samtools");
            Console.Out.WriteLine("samtools is not defined, default will be samtools. You can also modify it at file {0}.", Config.ConfigFilename);
            Config.Save();
          }

          if (null == this.BamFiles || 0 == this.BamFiles.Count)
          {
            ParsingErrors.Add("Bam files not defined.");
            return false;
          }

          if (2 != this.BamFiles.Count)
          {
            ParsingErrors.Add("Must be two indexed bam files for samtools mpileup.");
            return false;
          }

          foreach (var bamfile in this.BamFiles)
          {
            if (!File.Exists(bamfile))
            {
              ParsingErrors.Add(string.Format("Bam file is not exists {0}", bamfile));
              return false;
            }
          }
          if (null == this.GenomeFastaFile)
          {
            ParsingErrors.Add("Genome fasta file not defined.");
            return false;
          }
          if (!File.Exists(this.GenomeFastaFile))
          {
            ParsingErrors.Add(string.Format("Genome fasta file not exists {0}.", this.GenomeFastaFile));
            return false;
          }

          if (this.ThreadCount >= 2)
          {
            Console.WriteLine("Checking chromosome names for thread mode ...");
            if (this.ChromosomeNames == null || this.ChromosomeNames.Count == 0)
            {
              var fai = this.GenomeFastaFile + ".fai";
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
                this.ChromosomeNames = SequenceUtils.ReadFastaNames(this.GenomeFastaFile);
                if (this.ChromosomeNames.Count == 0)
                {
                  ParsingErrors.Add(string.Format("Genome fasta file doesn't contain chromosome names, {0}.", this.GenomeFastaFile));
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

          Console.Out.WriteLine("#mpileup bam files: " + this.BamFiles.Merge(","));
          Console.Out.WriteLine("#mpileup minimum read quality: " + this.MpileupMinimumReadQuality);
          Console.Out.WriteLine("#mpileup genome fasta: " + this.GenomeFastaFile);

          break;
        case DataSourceType.console:
          Console.Out.WriteLine("#mpileup from console.");
          break;
      }

      return true;
    }

    public string GetMpileupChromosomes()
    {
      if (this.ChromosomeNames != null && this.ChromosomeNames.Count > 0)
      {
        return this.ChromosomeNames.Merge(",");
      }
      else
      {
        return null;
      }
    }

    private bool PrepareOutputDirectory()
    {
      if (!Directory.Exists(this.OutputDirectory))
      {
        try
        {
          Directory.CreateDirectory(this.OutputDirectory);
        }
        catch (Exception ex)
        {
          ParsingErrors.Add(string.Format("Cannot create directory {0} : {1}", this.OutputDirectory, ex.Message));
          return false;
        }
      }

      if (!Directory.Exists(this.CandidatesDirectory))
      {
        try
        {
          Directory.CreateDirectory(this.CandidatesDirectory);
        }
        catch (Exception ex)
        {
          ParsingErrors.Add(string.Format("Cannot create directory {0} : {1}", this.CandidatesDirectory, ex.Message));
          return false;
        }
      }

      return true;
    }
  }
}
