using CQS.Genome.SomaticMutation;
using RCPA.Commandline;
using RCPA.Gui.Command;
using RCPA.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;

namespace CQS
{
  internal static class Program
  {
    private const int AttachParentProcess = -1;

    [DllImport("kernel32.dll")]
    private static extern bool AttachConsole(int dwProcessId);

    /// <summary>
    ///   The main entry point for the application.
    /// </summary>
    [STAThread]
    private static void Main(string[] args)
    {
      var commands = new ICommandLineCommand[]
      {
        new PipelineProcessorCommand(),
        //new FilterProcessorCommand(),
        new AnnotationProcessorCommand(),
        new ValidationProcessorCommand(),
        new ExtractProcessorCommand()
      }.ToDictionary(m => m.Name.ToLower());


      if (!SystemUtils.IsLinux)
      {
        AttachConsole(AttachParentProcess);
      }

      ICommandLineCommand command;
      if (args.Length == 0)
      {
        ShowUsage(commands);
      }
      else if (commands.TryGetValue(args[0].ToLower(), out command))
      {
        if (command.Process(args.Skip(1).ToArray()))
        {
          Console.WriteLine("Done!");
        }
        else
        {
          Console.Error.WriteLine("Failed!");
        }
      }
      else
      {
        Console.WriteLine("Error command " + args[0] + ".");
        ShowUsage(commands);
      }
    }

    private static void ShowUsage(Dictionary<string, ICommandLineCommand> commands)
    {
      Console.WriteLine(Constants.GetSqhVanderbiltTitle(GlmvcAssembly.Title, GlmvcAssembly.Version));
      Console.WriteLine("Those commands are available :");
      (from c in commands.Values
       select "\t" + c.Name + "\t" + c.Description).ToList().ForEach(Console.WriteLine);
    }
  }
}