using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CQS.Genome.Pileup;
using CQS.Statistics;
using System.IO;
using System.Diagnostics;
using CommandLine;
using RCPA.Utils;
using RCPA.Gui;
using System.Reflection;
using CQS.Genome.Statistics;
using System.Collections.Concurrent;
using System.Threading;

namespace wsmdetector
{
  class Program
  {
    static int totalCount = 0;
    static int minReadDepthFailed = 0;
    static int oneEventFailed = 0;
    static int minimumPercentageFailed = 0;
    static int groupFisherFailed = 0;
    static int positionFisherFailed = 0;
    static int strandFisherFailed = 0;
    static int candidateCount = 0;
    static int threadcount = 0;

    static List<BlockingCollection<string>> mpileuplinesList = new List<BlockingCollection<string>>();
    static BlockingCollection<string> mpileuplines = new BlockingCollection<string>();
    static ConcurrentQueue<FilterResult> filtered = new ConcurrentQueue<FilterResult>();
    static Options options = new Options();

    class FilterResult
    {
      public PileupItem Item { get; set; }
      public FisherExactTestResult Group { get; set; }
      public FisherExactTestResult Position { get; set; }
      public FisherExactTestResult Strand { get; set; }
    }

    class Processor
    {
      public int minReadDepthFailed = 0;
      public int oneEventFailed = 0;
      public int minimumPercentageFailed = 0;
      public int groupFisherFailed = 0;
      public int positionFisherFailed = 0;
      public int strandFisherFailed = 0;
      public int candidateCount = 0;

      PileupItemGroupTest groupTest = new PileupItemGroupTest();
      PileupItemPositionTest positionTest = new PileupItemPositionTest();
      PileupItemStrandTest strandTest = new PileupItemStrandTest();
      PileupItemPercentageTest percentageTest = new PileupItemPercentageTest(options.MinimumPercentageOfMinorAllele);
      PileupItemParser parser = GetParser();

      public FilterResult Parse(string line)
      {
        var item = parser.GetValue(line);
        if (item == null)
        {
          minReadDepthFailed++;
          return null;
        }

        //If the bases from all samples are same, ignore the entry.
        if (item.OnlyOneEvent())
        {
          oneEventFailed++;
          return null;
        }

        var fisherresult = item.InitializeTable();
        if (!percentageTest.Accept(fisherresult))
        {
          minimumPercentageFailed++;
          return null;
        }

        //group fisher exact test
        var groupFisher = groupTest.Test(fisherresult);
        if (groupFisher.PValue > options.PValue)
        {
          groupFisherFailed++;
          return null;
        }

        //position fisher exact test
        var positionFisher = new FisherExactTestResult();
        if (options.FilterPosition)
        {
          positionFisher = positionTest.Test(item);
          if (positionFisher.PValue <= options.PValue)
          {
            positionFisherFailed++;
            return null;
          }
        }

        //strand fisher exact test
        var strandFisher = new FisherExactTestResult();
        if (options.FilterStrand)
        {
          strandFisher = strandTest.Test(item);
          if (strandFisher.PValue <= options.PValue)
          {
            strandFisherFailed++;
            return null;
          }
        }

        candidateCount++;

        //get major and second alleles
        var bases = new HashSet<string>(new string[] { groupFisher.SucceedName, groupFisher.FailedName });

        //save to file
        PileupItemFile piFile = new PileupItemFile(bases);
        var filename = string.Format("{0}/{1}.wsm", options.OutputFileDirectory, piFile.GetFilename(item));
        piFile.WriteToFile(filename, item);

        return new FilterResult()
        {
          Item = item,
          Group = groupFisher,
          Position = positionFisher,
          Strand = strandFisher
        };
      }
    }

    static void FilterMpileupItem(Object stateInfo)
    {
      //var mpileuplines = stateInfo as BlockingCollection<string>;

      int mythreadindex = Interlocked.Increment(ref threadcount);
      Console.WriteLine("Start sub thread {0}", mythreadindex);
      try
      {
        Processor proc = new Processor();

        while (!mpileuplines.IsAddingCompleted)
        {
          Console.WriteLine("Waiting for data {0} ...", mythreadindex);
          Thread.Sleep(10);
          string line;
          while (mpileuplines.TryTake(out line))
          {
            var item = proc.Parse(line);

            if (item == null)
            {
              continue;
            }

            filtered.Enqueue(item);
          }
        }

        Interlocked.Add(ref minReadDepthFailed, proc.minReadDepthFailed);
        Interlocked.Add(ref oneEventFailed, proc.oneEventFailed);
        Interlocked.Add(ref minimumPercentageFailed, proc.minimumPercentageFailed);
        Interlocked.Add(ref groupFisherFailed, proc.groupFisherFailed);
        Interlocked.Add(ref positionFisherFailed, proc.positionFisherFailed);
        Interlocked.Add(ref strandFisherFailed, proc.strandFisherFailed);
        Interlocked.Add(ref candidateCount, proc.candidateCount);
      }
      finally
      {
        Interlocked.Decrement(ref threadcount);
        Console.WriteLine("Finished sub thread {0}", mythreadindex);
      }
    }

    private static PileupItemParser GetParser()
    {
      return new PileupItemParser(options.MinimumReadDepth, options.MinimumMappingQuality, options.IgnoreInsertionDeletion, true, options.IgnoreTerminalBase);
    }

    static void Main(string[] args)
    {
      ExecuteConfig config = ExecuteConfig.LoadExecuteConfig();
      if (config == null)
      {
        return;
      }

      Console.WriteLine("R:" + config.RExecute);
      Console.WriteLine("samtools:" + config.SamtoolsExecute);

      if (!options.ParseOptions(args, config))
      {
        return;
      }

      PileupItemParser parser = GetParser();
      PileupFile pfile = new PileupFile(parser);
      switch (options.From)
      {
        case SourceType.mpileup:
          pfile.Open(options.MpileupResultFile);
          break;
        case SourceType.bam:
          var proc = new Process
          {
            StartInfo = new ProcessStartInfo
            {
              FileName = config.SamtoolsExecute,
              Arguments = string.Format(" mpileup -q {0} -f {1} {2} {3} ", options.MpileupMinimumReadQuality, options.MpileupGenomeFile, options.MpileupBamFiles[0], options.MpileupBamFiles[1]),
              UseShellExecute = false,
              RedirectStandardOutput = true,
              CreateNoWindow = true
            }
          };

          Console.Out.WriteLine("running command : " + proc.StartInfo.FileName + " " + proc.StartInfo.Arguments);
          if (!proc.Start())
          {
            Console.Out.WriteLine("samtools mpileup cannot be started, check your parameters and ensure that samtools are available.");
            Console.Out.WriteLine(options.GetUsage());
            return;
          }

          pfile.Open(proc.StandardOutput);
          break;
        case SourceType.console:
          pfile.Open(Console.In);
          break;
      }

      Console.Out.WriteLine("#output directory: " + options.OutputDirectory);
      Console.Out.WriteLine("#minimum count: " + options.MinimumReadDepth.ToString());
      Console.Out.WriteLine("#minimum mapping quality: " + options.MinimumMappingQuality.ToString());
      Console.Out.WriteLine("#minimum percentage of minor allele: " + options.MinimumPercentageOfMinorAllele.ToString());
      Console.Out.WriteLine("#pvalue: " + options.PValue.ToString());
      Console.Out.WriteLine("#filter by position bias: " + options.FilterPosition.ToString());
      Console.Out.WriteLine("#filter by strand bias: " + options.FilterStrand.ToString());
      Console.Out.WriteLine("#thread count: " + options.ThreadCount.ToString());

      List<FilterResult> saved = new List<FilterResult>();

      Stopwatch watch = new Stopwatch();
      watch.Start();

      if (options.ThreadCount < 2)
      {
        using (pfile)
        {
          Processor proc = new Processor();

          string line;
          while ((line = pfile.ReadLine()) != null)
          //while ((item = pfile.Next("1", 48901870)) != null)
          {
            totalCount++;

            var item = proc.Parse(line);
            if (item == null)
            {
              continue;
            }

            saved.Add(item);
          }

          minReadDepthFailed = proc.minReadDepthFailed;
          oneEventFailed = proc.oneEventFailed;
          minimumPercentageFailed = proc.minimumPercentageFailed;
          groupFisherFailed = proc.groupFisherFailed;
          positionFisherFailed = proc.positionFisherFailed;
          strandFisherFailed = proc.strandFisherFailed;
          candidateCount = proc.candidateCount;
        }
      }
      else
      {
        var totalThread = options.ThreadCount - 1;
        for (int i = 0; i < totalThread; i++)
        {
          ThreadPool.QueueUserWorkItem(new WaitCallback(FilterMpileupItem));
          //var bc = new BlockingCollection<string>();
          //mpileuplinesList.Add(bc);
          //ThreadPool.QueueUserWorkItem(new WaitCallback(FilterMpileupItem), bc);
        }

        using (pfile)
        {
          string line;
          //int threadIndex = 0;
          while ((line = pfile.ReadLine()) != null)
          {
            totalCount++;

            mpileuplines.Add(line);

            //threadIndex = (threadIndex + 1) % totalThread;
            //mpileuplinesList[threadIndex].Add(line);
          }
          mpileuplines.CompleteAdding();
          //mpileuplinesList.ForEach(m => m.CompleteAdding());
        }

        while (threadcount > 0)
        {
          Thread.Sleep(100);
        }

        saved = filtered.ToList();
        saved.Sort((m1, m2) =>
        {
          var result = m1.Item.SequenceIdentifier.CompareTo(m2.Item.SequenceIdentifier);
          if (result == 0)
          {
            result = m1.Item.Position.CompareTo(m2.Item.Position);
          }
          return result;
        });
      }

      watch.Stop();
      Console.WriteLine("Time cost: {0}", watch.Elapsed);

      using (StreamWriter sw = new StreamWriter(options.OutputDirectory + "/candidates.summary"))
      {
        sw.WriteLine("chr\tloc\tref\tmajor_allele\tminor_allele\tsample_major_count\tsample_minor_count\tcontrol_major_count\tcontrol_minor_count\tfisher_group\tfisher_position\tfisher_strand");

        foreach (var res in saved)
        {
          sw.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9:0.0E00}\t{10:0.0E00}\t{11:0.0E00}",
            res.Item.SequenceIdentifier, res.Item.Position, res.Item.Nucleotide, res.Group.SucceedName, res.Group.FailedName,
            res.Group.SampleSucceedCount, res.Group.SampleFailedCount, res.Group.ControlSucceedCount, res.Group.ControlFailedCount, res.Group.PValue, res.Position.PValue, res.Position.PValue);
        }
      }

      using (StreamWriter sw = new StreamWriter(options.OutputDirectory + "/filter.summary"))
      {
        sw.WriteLine("Reason\tRemoved candidate\tRetained candidate");
        sw.WriteLine("total site\t\t{0}", totalCount);
        totalCount -= minReadDepthFailed;
        sw.WriteLine("minimum read depth failed\t{0}\t{1}", minReadDepthFailed, totalCount);
        totalCount -= oneEventFailed;
        sw.WriteLine("all same allele\t{0}\t{1}", oneEventFailed, totalCount);
        totalCount -= minimumPercentageFailed;
        sw.WriteLine("minimum percentage of minor allele failed\t{0}\t{1}", minimumPercentageFailed, totalCount);
        totalCount -= groupFisherFailed;
        sw.WriteLine("fisher exact test is not significant\t{0}\t{1}", groupFisherFailed, totalCount);
        if (options.FilterPosition)
        {
          totalCount -= positionFisherFailed;
          sw.WriteLine("fisher exact test of position is significant\t{0}\t{1}", positionFisherFailed, totalCount);
        }
        if (options.FilterStrand)
        {
          totalCount -= strandFisherFailed;
          sw.WriteLine("fisher exact test of strand is significant\t{0}\t{1}", strandFisherFailed, totalCount);
        }
      }

      var rproc = new Process
      {
        StartInfo = new ProcessStartInfo
        {
          //          FileName = @"C:\Program Files\R\R-2.15.2\bin\R.exe",
          FileName = config.RExecute,
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
      }
    }

  }
}
