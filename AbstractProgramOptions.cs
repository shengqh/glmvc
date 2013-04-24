using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CQS.Commandline;
using CQS;

namespace RSMC
{
  public abstract class AbstractProgramOptions : AbstractOptions
  {
    public CommandConfig Config { get; set; }
  }
}
