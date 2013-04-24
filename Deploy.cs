using System;
using NUnit.Framework;
using System.IO;
using RCPA.Utils;
using RCPA;

namespace Deployment
{
  [TestFixture]
  public class Deploy
  {
    [Test]
    public void DeployZip()
    {
      var appPath = Path.GetFullPath(".");
      Console.WriteLine(appPath);

      string dirName = "RSMC." + ThisAssembly.Version;
      string zipFile = DoZip(appPath, dirName);
    }

    private string DoZip(string appPath, string dirName)
    {
      string dirPath = new DirectoryInfo(appPath + "/" + dirName).FullName;

      CreateAndClearDirectory(dirPath);

      //string[] files = Directory.GetFiles(appPath);
      //foreach (var file in files)
      //{
      //  Console.WriteLine("\"{0}\",", Path.GetFileName(file));
      //}

      string[] filenames = new string[]{
        "CommandLine.dll",
        "CQS.Core.dll",
        "Meta.Numerics.dll",
        "RCPA.Core.dll",
        "RCPA.Proteomics.dll",
        "rsmc.exe",
        "rsmc.linux",
        "rsmc.r"
      };

      foreach (string filename in filenames)
      {
        File.Copy(appPath + "\\" + filename, dirPath + "\\" + filename);
      }

      string zipFile = dirName + ".zip";

      if (File.Exists(zipFile))
      {
        File.Delete(zipFile);
      }

      string zipCmdArgu = MyConvert.Format("a -tzip {0} {1}", zipFile, dirName);

      SystemUtils.Execute(@"C:\Program Files\7-Zip\7z.exe", zipCmdArgu, appPath);

      return zipFile;
    }

    private void CopyFile(string source, string target)
    {
      if (File.Exists(target))
      {
        File.Delete(target);
      }

      try
      {
        File.Copy(source, target);
      }
      catch (Exception)
      {
      }
    }

    private static void CopySubDirFiles(string sourcePath, string targetPath, string subDirName)
    {
      string scriptsPath = sourcePath + "/" + subDirName;
      string scriptsDir = targetPath + "/" + subDirName;

      if (Directory.Exists(scriptsPath))
      {
        CreateAndClearDirectory(scriptsDir);

        string[] scriptFiles = Directory.GetFiles(scriptsPath);
        foreach (string scriptFile in scriptFiles)
        {
          File.Copy(scriptFile, scriptsDir + "\\" + new FileInfo(scriptFile).Name);
        }
      }
    }

    private static void CreateAndClearDirectory(string dirPath)
    {
      if (Directory.Exists(dirPath))
      {
        string[] files = Directory.GetFiles(dirPath);
        foreach (string file in files)
        {
          File.Delete(file);
        }
      }
      else
      {
        Directory.CreateDirectory(dirPath);
      }
    }
  }
}
