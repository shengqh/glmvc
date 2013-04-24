using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CommandLine;
using CommandLine.Text;

namespace RSMC
{
  public class Options
  {
    [VerbOption("pileup", HelpText = "Initialize candidates from samtools mpileup result.")]
    public PileupOptions InitVerb { get; set; }

    [VerbOption("filter", HelpText = "Filter candidates by logistic regression model.")]
    public FilterOptions FilterVerb { get; set; }

    [VerbOption("annotation", HelpText = "Annotate mutation using varies tools.")]
    public AnnotationOptions AnnotationVerb { get; set; }

    [VerbOption("all", HelpText = "pileup/filter/annotate data")]
    public AllOptions AllVerb { get; set; }

    [HelpVerbOption]
    public string GetUsage(string verb)
    {
      return HelpText.AutoBuild(this, verb);
    }
  }
}
