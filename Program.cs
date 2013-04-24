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
using RCPA.Seq;
using CQS.Genome.Annotation;
using CQS;

namespace RSMC
{
  class Program
  {
    static void Main(string[] args)
    {
      var config = new CommandConfig();
      config.Load();

      if (args.Length < 1)
      {
        Console.WriteLine(new Options().GetUsage(null));
        Environment.Exit(-1);
      }

      string invokedVerb = null;
      AbstractProgramOptions invokedVerbInstance = null;
      Options options = null;
      try
      {
        options = CommandLine.Parser.Default.ParseArguments<Options>(args,
          (verb, subOptions) =>
          {
            if (subOptions != null)
            {
              invokedVerb = verb;
              invokedVerbInstance = (AbstractProgramOptions)subOptions;
              invokedVerbInstance.Config = config;
              invokedVerbInstance.IsPileup = false;
            }
          },
          () =>
          {
            Environment.Exit(-1);
          }
        );
      }
      catch (Exception ex)
      {
        Console.Out.WriteLine(ex.Message);
        Environment.Exit(-1);
      }

      if (invokedVerb == "init")
      {
        var initOptions = (PileupOptions)invokedVerbInstance;
        if (!initOptions.PrepareOptions())
        {
          Console.Out.WriteLine(initOptions.GetUsage());
          Environment.Exit(-1);
        }

        var samtools = new PileupProcessor(initOptions);
        if (!samtools.Process())
        {
          Environment.Exit(-1);
        }
      }
      else if (invokedVerb == "filter")
      {
        var filterOptions = (FilterOptions)invokedVerbInstance;
        if (!filterOptions.PrepareOptions())
        {
          Console.Out.WriteLine(filterOptions.GetUsage());
          Environment.Exit(-1);
        }

        var rtools = new FilterProcessor(filterOptions);
        if (!rtools.Process())
        {
          Environment.Exit(-1);
        }
      }
      else if (invokedVerb == "annotation")
      {
        var annotationOptions = (AnnotationOptions)invokedVerbInstance;
        if (!annotationOptions.PrepareOptions())
        {
          Console.Out.WriteLine(annotationOptions.GetUsage());
          Environment.Exit(-1);
        }

        var annovar = new AnnotationProcessor(annotationOptions);
        if (!annovar.Process())
        {
          Environment.Exit(-1);
        }
      }
      else if (invokedVerb == "all")
      {
        var allOptions = (AllOptions)invokedVerbInstance;

        //check paramater first
        if (!allOptions.PrepareOptions())
        {
          Console.Out.WriteLine(allOptions.GetUsage());
          Environment.Exit(-1);
        }
        FilterOptions filterOptions = allOptions.GetFilterOptions();
        if (!filterOptions.PrepareOptions())
        {
          allOptions.ParsingErrors.AddRange(filterOptions.ParsingErrors);
          Console.Out.WriteLine(allOptions.GetUsage());
          Environment.Exit(-1);
        }
        AnnotationOptions annotationOptions = allOptions.GetAnnotationOptions();
        if (!annotationOptions.PrepareOptions())
        {
          allOptions.ParsingErrors.AddRange(annotationOptions.ParsingErrors);
          Console.Out.WriteLine(allOptions.GetUsage());
          Environment.Exit(-1);
        }

        //run initialize candidates
        var samtools = new PileupProcessor(allOptions);
        if (!samtools.Process())
        {
          Environment.Exit(-1);
        }

        //check the result exists
        filterOptions.IsPileup = false;
        if (!filterOptions.PrepareOptions())
        {
          foreach (var line in filterOptions.ParsingErrors)
          {
            Console.WriteLine("errors : " + line);
          }
          Environment.Exit(-1);
        }

        var filter = new FilterProcessor(filterOptions);
        if (!filter.Process())
        {
          Environment.Exit(-1);
        }

        annotationOptions.IsPileup = false;
        if (!annotationOptions.PrepareOptions())
        {
          foreach (var line in annotationOptions.ParsingErrors)
          {
            Console.WriteLine("errors : " + line);
          }
          Environment.Exit(-1);
        }

        var annotation = new AnnotationProcessor(annotationOptions);
        if (!annotation.Process())
        {
          Environment.Exit(-1);
        }
      }
      else
      {
        Console.WriteLine("I don't know the command " + invokedVerb);
        Environment.Exit(-1);
      }

      Console.WriteLine("All done!");
    }
  }
}