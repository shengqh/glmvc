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
    [VerbOption("init", HelpText = "Initialize candidates from samtools mpileup result.")]
    public InitOptions InitVerb { get; set; }

    [VerbOption("filter", HelpText = "Filter candidates by logistic regression model.")]
    public FilterOptions FilterVerb { get; set; }

    [VerbOption("annotation", HelpText = "Annotate mutation using varies tools.")]
    public AnnotationOptions AnnotationVerb { get; set; }

    [VerbOption("all", HelpText = "Parse/filter/annotate data")]
    public AllOptions AllVerb { get; set; }

    [HelpVerbOption]
    public string GetUsage(string verb)
    {
      return HelpText.AutoBuild(this, verb);
    }
  }
}
