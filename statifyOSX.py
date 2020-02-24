import os, sys, string, subprocess

libsRequired = set()
updateFiles = False
processLibs = False
OSTYPE = "Linux"
MACHINETYPE="x86_64"
OSVERSION = "1"
GCC_VERSION="0"

def getGCCVersion():
   global GCC_VERSION

   try:
      GCC_VERSION=(getCommandOutput("gcc --version|grep gcc|awk '{print $NF}' |awk -F \".\" '{print $1$2}'", False))
   except:
      try:
         GCC_VERSION=(getCommandOutput("gcc --version|grep gcc|awk '{print $3}' |awk -F \".\" '{print $1$2}'", False))
      except:
         print ("Warning: cannot determine GCC version")
   GCC_VERSION = "gcc" + GCC_VERSION.decode("utf-8")

def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1")

def getCommandOutput(theCommand, checkForStderr):
    p = subprocess.Popen(theCommand, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (checkStdout, checkStderr) = p.communicate()
    if checkForStderr and checkStderr != "":
       return ""
    else:
       return checkStdout.strip()

def getMachineType():
   global OSTYPE
   global OSVERSION
   global MACHINETYPE
   getGCCVersion()

   p = subprocess.Popen("echo `uname`", shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   (checkStdout, checkStderr) = p.communicate()
   if checkStderr.decode("utf-8") != "":
      print ("Warning: Cannot determine OS, defaulting to %s"%(OSTYPE))
   else:
      OSTYPE = checkStdout.decode("utf-8").strip()

   p = subprocess.Popen("echo `uname -r`", shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   (checkStdout, checkStderr) = p.communicate()
   if checkStderr.decode("utf-8") != "":
      print ("Warning: Cannot determine OS version, defaulting to %s"%(OSVERSION))
   else:
      OSVERSION = checkStdout.decode("utf-8").strip()

   p = subprocess.Popen("echo `uname -m`", shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   (checkStdout, checkStderr) = p.communicate()
   if checkStderr.decode("utf-8") != "":
      print ("Warning: Cannot determine system type, defaulting to %s"%(MACHINETYPE))
   else:
      MACHINETYPE = checkStdout.decode("utf-8").strip()

def is_binary(filename):
   fin = open(filename, 'rb')
   try:
       CHUNKSIZE = 1024
       while 1:
           chunk = fin.read(CHUNKSIZE)
           if b'\0' in chunk: # found null byte
               return True
           if len(chunk) < CHUNKSIZE:
               break # done
   finally:
       fin.close()

   return False

def statifyBin(binName, relativePath):
   global updateFiles
   global libsRequired

   print ("Processing bin %s with relative path %s"%(binName, relativePath))

   cmdOut = ""
   if OSTYPE == "Darwin":
      cmdOut = getCommandOutput("otool -L %s"%(binName), False)
   elif OSTYPE == "Linux":
      cmdOut = getCommandOutput("ldd %s"%(binName), False)
   else:
      print ("Unknown OS TYPE")
      return

   if binName.endswith("dylib") and not processLibs:
      return

   binDir = os.path.dirname(binName)
   libDir = os.path.abspath(binDir + os.sep + relativePath)
   isNotStatic = False
   print ("For bin %s in lib %s processing %s"%(binName, libDir, cmdOut))
   for line in cmdOut.decode("utf-8").split("\n"):
      if "compatibility" in line:
         if "@executable_path" in line:
            fileToCheck = os.path.abspath(binDir + os.sep + line.replace("@executable_path", "").split()[0].strip())
            if not os.path.exists(fileToCheck):
               print ("Error file %s doesnt exist for %s"%(fileToCheck, binName))
               sys.exit(1)
            continue
         if "System/Library/" in line:
            continue
         if "libSystem.B.dylib" in line:
            continue
         if "/usr/lib/" in line and "libmpi" not in line and "libopen" not in line:
            continue
         if "/usr/local/lib" in line and "lbmpi" not in line and "libopen" not in line:
            continue
         libName = line.strip().split()[0]
         libFile = os.path.basename(libName)
         if "/gcc" in libName and GCC_VERSION not in libName:
            print ("Warning: skipped adding dependency for %s (%s) due to version conflict with current gcc %s"%(binName, libName, GCC_VERSION))
         else:
            libsRequired.add(os.path.abspath(libName))
         if updateFiles and OSTYPE == "Darwin":
            print ("Running command install_name_tool -change %s @executable_path%s%s%s%s %s"%(libName, os.sep, relativePath,os.sep, libFile, binName))
            os.system("install_name_tool -change %s @executable_path%s%s%s%s %s"%(libName, os.sep, relativePath,os.sep, libFile, binName))
         isNotStatic = True
      elif "=>" in line and "(0x0" in line and libDir not in line:
         libIDS = line.split()
         libName = libIDS[len(libIDS)-2]
         if "(0x)" in libName or "=>" in libName:
            continue
         libsRequired.add(os.path.abspath(libName))

   if isNotStatic:
       print ("File %s is not statically compiled!"%(binName))

# traverse root directory, and list directories as dirs and files as files
def walkdir(dirname, libDir):
    for cur, dirs, files in os.walk(dirname):
       relativePath=""
       #if (dirname != libDir):
       for i in range(0,len(cur.split(os.sep))):
          relativePath = "../%s"%(relativePath)
       #else:
       #   relativePath=""

       for f in files:
          fullPath = cur + os.sep + f
          if os.path.islink(fullPath) or not os.path.exists(fullPath):
             continue
          if is_binary(fullPath):
             statifyBin(fullPath, relativePath + os.sep + libDir)
    return

def main():
   global updateFiles
   global processLibs
   global libsRequired
   libs = "lib"
   dirName = "."

   getMachineType()
   if (len(sys.argv) > 4):
      processLibs = str2bool(sys.argv[4])
   if (len(sys.argv) > 3):
      updateFiles = str2bool(sys.argv[3])
   if (len(sys.argv) > 2):
      libs = sys.argv[2]
   if (len(sys.argv) > 1):
      dirName = os.path.relpath(sys.argv[1])

   os.system("mkdir -p %s"%(libs))
   walkdir(dirName, libs)
   print ("Dependencies found:\n%s"%(libsRequired))

   for lib in libsRequired:
      fileName = os.path.basename(lib)
      if not os.path.exists("%s%s%s"%(libs, os.sep, fileName)):
         os.system("cp %s %s/"%(lib, libs))
      else:
         print ("Warning: did not copy %s because same file name already exists"%(lib))

if __name__ == '__main__':
    main()
