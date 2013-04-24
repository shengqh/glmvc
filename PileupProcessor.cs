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
  public class PileupProcessor
  {
    private PileupOptions options;

    public PileupProcessor(PileupOptions options)
    {
      this.options = options;
    }

    private int totalCount = 0;
    private int minReadDepthFailed = 0;
    private int oneEventFailed = 0;
    private int minimumPercentageFailed = 0;
    private int groupFisherFailed = 0;
    private int positionFisherFailed = 0;
    private int strandFisherFailed = 0;
    private int candidateCount = 0;
    private int threadcount = 0;

    private int startedThreadCount = 0;

    private bool? samtools_ok = null;

    private int thread_failed = 0;

    private ConcurrentQueue<FilterResult> filtered = new ConcurrentQueue<FilterResult>();

    private class FilterResult
    {
      public PileupItem Item { get; set; }
      public FisherExactTestResult Group { get; set; }
      public FisherExactTestResult Position { get; set; }
      public FisherExactTestResult Strand { get; set; }
    }

    private class Processor
    {
      public int minReadDepthFailed = 0;
      public int oneEventFailed = 0;
      public int minimumPercentageFailed = 0;
      public int groupFisherFailed = 0;
      public int positionFisherFailed = 0;
      public int strandFisherFailed = 0;
      public int candidateCount = 0;

      private PileupOptions options;

      public Processor(PileupOptions options)
      {
        this.options = options;
        this.rdFilter = new PileupItemReadDepthFilter(options.MinimumReadDepth, options.MinimumBaseQuality);
        this.groupTest = new PileupItemGroupTest();
        this.positionTest = new PileupItemPositionTest();
        this.strandTest = new PileupItemStrandTest();
        this.percentageTest = new PileupItemPercentageTest(options.MinimumPercentageOfMinorAllele);
        this.parser = options.GetPileupItemParser();
        this.scoreFilter = m => m.Score >= this.options.MinimumBaseQuality;
      }

      private PileupItemReadDepthFilter rdFilter;
      private PileupItemGroupTest groupTest;
      private PileupItemPositionTest positionTest;
      private PileupItemStrandTest strandTest;
      private PileupItemPercentageTest percentageTest;
      private PileupItemParser parser;
      private Func<PileupBase, bool> scoreFilter;

      public FilterResult Parse(string line)
      {
        var item = parser.GetValue(line);
        if (item == null)
        {
          minReadDepthFailed++;
          return null;
        }

        if (!rdFilter.Accept(item))
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

        //only sample2 (tumor sample) was tested.
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
        var filename = string.Format("{0}/{1}.wsm", options.CandidatesDirectory, piFile.GetFilename(item));
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

    private void ParallelChromosome(Object stateInfo)
    {
      ConcurrentQueue<string> chromosomes = stateInfo as ConcurrentQueue<string>;

      Interlocked.Increment(ref startedThreadCount);

      int mythreadindex = Interlocked.Increment(ref threadcount);
      Console.WriteLine("Start sub thread {0}", mythreadindex);
      try
      {
        Processor proc = new Processor(this.options);
        PileupItemParser parser = options.GetPileupItemParser();

        int localTotalCount = 0;

        while (!chromosomes.IsEmpty)
        {
          string chromosomeName;
          if (!chromosomes.TryDequeue(out chromosomeName))
          {
            Thread.Sleep(10);
            continue;
          }

          Console.WriteLine("Processing chromosome {0} in thread {1}", chromosomeName, mythreadindex);
          var process = ExecuteSamtools(options.GetSamtoolsCommand(), chromosomeName);
          if (process == null)
          {
            return;
          }

          try
          {
            using (PileupFile pfile = new PileupFile(parser))
            {
              pfile.Open(process.StandardOutput);
              string line;
              while ((line = pfile.ReadLine()) != null)
              {
                localTotalCount++;

                if (thread_failed > 0)
                {
                  process.Kill();
                  return;
                }

                try
                {
                  var item = proc.Parse(line);
                  if (item == null)
                  {
                    continue;
                  }

                  filtered.Enqueue(item);
                }
                catch (Exception ex)
                {
                  Console.WriteLine("Parsing mpileup result error : {0}\n{1}", ex.Message, line);
                  Interlocked.Increment(ref thread_failed);
                  return;
                }
              }
            }
          }
          finally
          {
            try { process.Kill(); }
            catch (Exception) { }
          }
        }

        Interlocked.Add(ref totalCount, localTotalCount);
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

    private void ParallelMpileupItem(Object stateInfo)
    {
      var mpileuplines = stateInfo as BlockingCollection<string>;

      int mythreadindex = Interlocked.Increment(ref threadcount);
      Console.WriteLine("Start sub thread {0} at {1}", mythreadindex, DateTime.Now);
      try
      {
        Processor proc = new Processor(this.options);

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
        Console.WriteLine("Finished sub thread {0} at {1}", mythreadindex, DateTime.Now);
      }
    }

    public bool Process()
    {
      Console.Out.WriteLine("initialize process started at {0}", DateTime.Now);

      Console.Out.WriteLine("#output directory: " + options.OutputDirectory);
      Console.Out.WriteLine("#minimum count: " + options.MinimumReadDepth.ToString());
      Console.Out.WriteLine("#minimum mapping quality: " + options.MinimumBaseQuality.ToString());
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
        Console.WriteLine("Single thread mode ...");
        PileupItemParser parser = options.GetPileupItemParser();
        PileupFile pfile = new PileupFile(parser);
        switch (options.From)
        {
          case PileupOptions.DataSourceType.mpileup:
            pfile.Open(options.MpileupFile);
            break;
          case PileupOptions.DataSourceType.bam:
            Process proc = ExecuteSamtools(options.GetSamtoolsCommand(), options.GetMpileupChromosomes());
            if (proc == null)
            {
              return false;
            }

            pfile.Open(proc.StandardOutput);
            pfile.Samtools = proc;
            break;
          case PileupOptions.DataSourceType.console:
            pfile.Open(Console.In);
            break;
        }

        using (pfile)
        {
          try
          {
            Processor proc = new Processor(this.options);

            string line;
            while ((line = pfile.ReadLine()) != null)
            //while ((item = pfile.Next("1", 48901870)) != null)
            {
              totalCount++;

              try
              {
                var item = proc.Parse(line);
                if (item == null)
                {
                  continue;
                }

                saved.Add(item);
              }
              catch (Exception ex)
              {
                Console.WriteLine("parsing error {0}\n{1}", ex.Message, line);
                return false;
              }
            }

            minReadDepthFailed = proc.minReadDepthFailed;
            oneEventFailed = proc.oneEventFailed;
            minimumPercentageFailed = proc.minimumPercentageFailed;
            groupFisherFailed = proc.groupFisherFailed;
            positionFisherFailed = proc.positionFisherFailed;
            strandFisherFailed = proc.strandFisherFailed;
            candidateCount = proc.candidateCount;
          }
          finally
          {
            if (pfile.Samtools != null)
            {
              try { pfile.Samtools.Kill(); }
              catch (Exception) { }
            }
          }
        }
      }
      else
      {
        var totalThread = options.ThreadCount - 1;
        if (options.From == PileupOptions.DataSourceType.bam)
        {
          Console.WriteLine("Multiple thread mode, parallel by chromosome ...");
          ConcurrentQueue<string> chromosomes = new ConcurrentQueue<string>();
          foreach (var chr in options.ChromosomeNames)
          {
            chromosomes.Enqueue(chr);
          }

          samtools_ok = null;
          ThreadPool.QueueUserWorkItem(new WaitCallback(ParallelChromosome), chromosomes);
          while (samtools_ok == null)
          {
            Thread.Sleep(100);
          }

          if (samtools_ok == false)
          {
            return false;
          }

          for (int i = 0; i < totalThread - 1; i++)
          {
            ThreadPool.QueueUserWorkItem(new WaitCallback(ParallelChromosome), chromosomes);
          }

          while (startedThreadCount == 0)
          {
            Thread.Sleep(100);
          }
        }
        else
        {
          Console.WriteLine("Multiple thread mode, parallel by mpileup item ...");
          BlockingCollection<string> mpileuplines = new BlockingCollection<string>();
          for (int i = 0; i < totalThread; i++)
          {
            ThreadPool.QueueUserWorkItem(new WaitCallback(ParallelMpileupItem), mpileuplines);
          }

          PileupItemParser parser = options.GetPileupItemParser();
          using (PileupFile pfile = new PileupFile(parser))
          {
            if (options.From == PileupOptions.DataSourceType.mpileup)
            {
              pfile.Open(options.MpileupFile);
            }
            else
            {
              pfile.Open(Console.In);
            }

            string line;
            while ((line = pfile.ReadLine()) != null)
            {
              totalCount++;

              mpileuplines.Add(line);
            }
            mpileuplines.CompleteAdding();
          }
        }

        while (threadcount > 0)
        {
          Thread.Sleep(100);
        }

        if (thread_failed > 0)
        {
          return false;
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
      Console.Out.WriteLine("initialize process ended at {0}, cost {1}", DateTime.Now, watch.Elapsed);

      using (StreamWriter sw = new StreamWriter(options.OutputDirectory + "/candidates.summary"))
      {
        sw.WriteLine("chr\tloc\tref\tmajor_allele\tminor_allele\tsample_major_count\tsample_minor_count\tcontrol_major_count\tcontrol_minor_count\tfisher_group\tfisher_position\tfisher_strand");

        foreach (var res in saved)
        {
          sw.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9:0.0E00}\t{10:0.0E00}\t{11:0.0E00}",
            res.Item.SequenceIdentifier, res.Item.Position, res.Item.Nucleotide, res.Group.SucceedName, res.Group.FailedName,
            res.Group.SucceedCount1, res.Group.FailedCount1, res.Group.SucceedCount2, res.Group.FailedCount2, res.Group.PValue, res.Position.PValue, res.Position.PValue);
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

      return true;
    }

    private Process ExecuteSamtools(string samtoolsCommand, string chromosome)
    {
      var chr = string.IsNullOrEmpty(chromosome) ? "" : " -r " + chromosome;
      var result = new Process
      {
        StartInfo = new ProcessStartInfo
        {
          FileName = samtoolsCommand,
          Arguments = string.Format(" mpileup -q {0} {1} -f {2} {3} {4} ", options.MpileupMinimumReadQuality, chr, options.GenomeFastaFile, options.BamFiles[0], options.BamFiles[1]),
          UseShellExecute = false,
          RedirectStandardOutput = true,
          CreateNoWindow = true
        }
      };

      Console.Out.WriteLine("running command : " + result.StartInfo.FileName + " " + result.StartInfo.Arguments);
      try
      {
        if (!result.Start())
        {
          Console.Out.WriteLine("samtools mpileup cannot be started, check your parameters and ensure that samtools are available.");
          samtools_ok = false;
          return null;
        }
      }
      catch (Exception ex)
      {
        Console.Out.WriteLine("samtools mpileup cannot be started, check your parameters and ensure that samtools are available : {0}", ex.Message);
        samtools_ok = false;
        return null;
      }

      samtools_ok = true;
      return result;
    }
  }
}