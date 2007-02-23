#!/usr/bin/env ruby

######################################################################
# test_ca.rb
# 
# $Id: test_ca.rb,v 1.4 2007-02-23 21:27:59 catmandew Exp $
#
# Terms:
#   test spec - combination of version, run script, config file, recipe script,
#               and any other parameters used to run a version
#
# Inputs:
#   1. set of test specs
#        includes
#          version path
#          path/name of run script
#          flag and path/name of run script's config file
#          flag and path/name of run script's recipe script (if relevant)
#          additional flags
#      designation of one as baseline for comparison of results
#   2. set of genomes to run on
#   3. path in which to put assembler outputs
#   4. path of qcCompare and tampaCompare for inter-assembly comparisons
#
# Outputs:
#   For each genome,
#     1. location, run time, success/failure of each test spec
#     2. for each failure, specific failure point (if it can be determined)
#     3. for each success:
#          run time
#          qc comparison relative to baseline
#          tampa comparison relative to baseline
#          judgment as to whether superior or inferior relative to baseline(?)
#          
#
######################################################################

require 'getoptlong'

STUBBING = true
# STUBBING = false

TESTING = true
# TESTING = false

Dependencies = {:pullfrag => "/local/common/pullfrag",
                :html2text => "/usr/bin/html2text"}

# files needed for carun, created by pullfrag
InputFileSuffixes = {:frag => "frg",
  :insert => "insert",
  :project => "project",
  :seq_features => "seq.features"}

TampaFileSuffixes = {:intra_breakpoints => "intra.breakpoints.tampa",
  :intra_summary => "intra.summary.tampa",
  :inter_breakpoints => "inter.breakpoints.tampa",
  :inter_summary => "inter.summary.tampa"}

AServerConsole = "http://assemblyconsole.tigr.org:8080/AserverConsole"
AServerRootDir = "/local/aserver"
RequestString = "requestInfo?id="
CheckString = "Current State: <INITIALIZING|RUNNING|ERROR|FINISHED>"

######################################################################
# print_help
def print_help(message = nil)
  warn "\n!!! #{message} !!!" if(message)

  spaces = ""
  basename = File.basename($0)
  basename.length.times{spaces += " "}

  printf(STDERR, "
Test different versions of the Celera Assembler on different genomes
and compare results.

Usage: #{basename}
           [--help|-h]
           [--verbose|-v]
           <--testspec|-t testSpecFile1> ...
           <--Genomes|-G genomesFile1>  ...
           <--genome|-g genome1>  ...
           [-comp|-c comparisonScriptsPath]

  testSpecFiles:
    At least one must be specified. These are tab-delimited text files
    that specify one test specification per line. Each specification
    is defined by the following fields.

      field 1: CA version name (e.g., 3.20)
      field 2: label for test spec (e.g., withOBT)
        The concatenation of CA version name and label with white space
        removed must be unique over the set of all test specifications
        for this run of #{$0}.
      field 3: full path to assembler binaries
               (e.g., /local/devel/SE/work/me/wgs-assembler/Linux64/bin)
      field 4: full path and name of run script
               name prefix must either carun or run_CA
      field 5: full path and name of config file
      field 6: full path and name of recipe script (blank if not carun)
      field 7: any extra parameters to be passed to run script
        (There is no need to add -p CA_BIN=<binaries path> for carun)

    For each genome, the earliest test specification listed that runs
    to completion will be treated as the 'baseline' version for purposes
    of results comparison.

  genomeFiles and genomes:
    A genomesFile lists one genome per line.
    Either at least one genomeFile must be specified listing at least one
    genome, or at least one genomes must be specified on the command line.
    All genomes are identified by their database identifier. For example,
    the genome identifier for D. mojavensis Wolbachia is gdmo.
    (Use the getdb program for assistance.)

  comparisonScriptsPath is the full path of the directory in which
    the qcCompare and tampaCompare scripts reside. Generally, these
    will be with the most recent version of assembler being tested


Processing:
  A subdirectory with name ca_test_<date>_<ttt> will be created to hold
  various files including a log file and results files.

  For each genome, #{$0} does the following:
    1. Creates a subdirectory to hold assembler input files, if needed,
    2. Runs pullfrag to create assembler input files, if needed,
    3. Runs each test on the grid,
    4. Runs TAMPA (if necessary), qcCompare, and tampaCompare to compare
       test results, and
    5. Creates a <genome> results file.

Outputs:
  Log output is written to log.txt, which lists the following:
    1. Status of test script steps,
    2. Success or failure of each test spec/genome combination

  A report is written for each genome to <genome>_results.txt,
  which lists the following for each successful test spec:
    1. qc comparison relative to the baseline test spec
       (see above for discussion of baseline), and
    2. TAMPA comparison relative to baseline.

")

  exit 1
end
# end of print_help function
######################################################################


######################################################################
# augment the File class with a method to get to a real dir or path
# through however many links may be involved
class File
  def File.original(name)
    return nil if(!File.exists?(name))
    true_path = name
    true_path = File.readlink(true_path) while(File.ftype(true_path) == "link")
    File.exists?(true_path) ? File.expand_path(true_path) : nil
  end  

  def File.get_name_and_original(full_path)
    return [File.split(full_path)[-1], File.original(full_path)]
  end
end # end of File class modification
######################################################################


######################################################################
# augment the Dir class with a method to get to create a subdir
# or verify it exists and is a directory
class Dir
  def Dir.create_sub_dir(name)
    Dir.mkdir(name) if(!File.original(name))
      
    if(!File.original(name) || File.ftype(File.original(name)) != "directory")
      return nil
    else
      return File.original(name)
    end
  end

  def Dir.create_path(name)
    if(!File.original(name))
      dirs = name.split("/")
      path = ""
      dirs.each do |dir|
        path = path + "/" + dir
        return nil if(!Dir.create_sub_dir(path))
      end
    elsif(File.ftype(File.original(name)) != "directory")
      return nil
    end
    return File.original(name)
  end
end # end of Dir class modification
######################################################################


######################################################################
# Logging class
class Logger
  def initialize(filename, be_verbose = false)
    @filename = filename
    @fo = File.open(filename, "w")
    @be_verbose = be_verbose
  end

  def log(string, time_stamp = false)
    time = (time_stamp ? "#{Time.now.strftime("%Y%m%d %H:%M:%S")}\t" : "")
    printf(@fo, "#{time}#{string}\n")
    @fo.flush
    printf(STDERR, "#{time}#{string}\n") if(@be_verbose)
  end
end # end of log file class
######################################################################


######################################################################
# Assembler Test Spec Class
#
# Spec file format is tab-delimited text, one spec per line. Fields are:
# 1. assembler version label
# 2. spec label
# 3. assembler binary path
# 4. script/config file for carun or run_CA
# 5. additional parameters to pass to carun or run_CA
#
class AssemblerTestSpec
  attr_reader :version, :label
  attr_reader :binary_path
  attr_reader :script, :script_type
  attr_reader :config_file
  attr_reader :recipe_script
  attr_reader :params
  attr_reader :line

  def initialize(spec_line)
    # initialize to bad values
    @line = spec_line
    @version = @label = ""
    @binary_path = @script = @config_file = @recipe_script = @params = nil
    @script_type = "invalid"

    # split the line on tabs, clean up, check number of fields
    fs = spec_line.chomp.split(/\t/,-1).map{|val| val.strip}

    return if(fs.size != 7)
    # set members
    @version = fs[0]
    @label = fs[1]
    @binary_path = fs[2]
    @script = fs[3]
    @config_file = fs[4]
    @recipe_script = fs[5]
    @params = fs[6]

    # check that script is either carun or run_CA
    script_basename = File.basename(@script)
    if(script_basename =~ /^carun/)
      @script_type = "carun"
    elsif(script_basename =~ /^run_CA/)
      @script_type = "run_CA"
      @recipe_script = "N/A"
    end
  end

  def version_label
    "#{@version}_#{@label}"
  end

  def key
    "Binaries: #{@binary_path}, Script: #{@script}, Recipe script: #{@recipe_script}, Params: #{@params}"
  end

  def errors
    errors = []
    errors << "No version" if(@version == "")
    errors << "No label" if(@label == "")
    errors << "No binary path" if(@binary_path == nil)
    errors << "Binary path does not exist" if(!File.exists?(@binary_path) ||
                                              File.ftype(@binary_path) != "directory")
    errors << "No script" if(@script == nil)
    errors << "Script does not exist" if(!File.exists?(@script))
    errors << "No parameters" if(@params == nil)
    errors << "No script type" if(@script_type == "invalid")
    errors << "No recipe script" if(@recipe_script == "" &&
                                    @script_type == "carun")
    errors << "Recipe script does not exist" if(@script_type == "carun" &&
                                                !File.exists?(@recipe_script))
                                                
    errors << "Whitespace in version" if(@version =~ /\s/)
    errors << "Whitespace in label" if(@label =~ /\s/)
    return errors.join("\n") if(errors.size > 0)
    return nil
  end
end # end of assembler test spec class
######################################################################


######################################################################
# Assembler Test Class
class AssemblerTest
  
  attr_reader :test_spec, :genome, :test_dir
  attr_reader :command, :config_file
  attr_reader :request_file
  attr_reader :aserver_dir

  def initialize(test_spec, genome, test_dir)
    @test_spec = test_spec
    @genome = genome
    @test_dir = test_dir
    @status = "NOT READY"
    @config_file = @test_spec.config_file
    @command = nil
    @request_file = "#{@test_dir}/#{test_spec.version_label}_#{@genome.name}_request.txt"
    @request_id = nil
    @aserver_dir = nil
  end

  def set_up
    return if(@status != "NOT READY")

    # for run_CA, need to modify config file
    case @test_spec.script_type
    when "run_CA"
      @config_file = "#{@test_dir}/run_CA_#{@test_spec.version_label}.config"
      if(!File.exists?(@config_file))
        root_dir = File.dirname(@test_spec.binary_path)
        ofp = File.open(@config_file, "w")
        open(@test_spec.config_file).each do |line|
          line.chomp!
          line = "root = #{root_dir}" if(line =~ /^\s*root\s*=/)
          printf(ofp, "#{line}\n")
        end # end of reading config file
        ofp.close
      end # end of check for test spec config file
      @command = "#{@test_spec.script} -C #{@config_file} -nowait -noedit -noinsert -nojoin -clean 1 #{@genome.dir}/#{@genome.name}.frg > #{@request_file}"
    when "carun"
      @command = "#{@test_spec.script} -config #{@config_file} -noCopy -R #{@test_spec.recipe_script} -p CA_BIN=#{@test_spec.binary_path} #{@genome.dir}/#{@genome.name}.frg > #{@request_file}"
    end # end of case check for carun/run_CA
    @status = "READY"
  end # end of set_up method

  def run
    return false if(@status != "READY")
    system(@command)
    if($? != 0)
      @status = "INITIATION ERROR"
      return false
    else
      @status = "INITIALIZING"
      return true
    end # end of if command ran ok
  end # end of run method

  def request_id
    return @request_id if(@request_id)
    return nil if(!@request_file || !File.exists?(@request_file))
    pattern = /The request id is ([0-9]+)/
    open(@request_file).each do |line|
      if(line =~ pattern)
        @request_id = line.match(pattern).to_a[1]
        @aserver_dir = AServerRootDir + "/" + @request_id
        break
      end # end of check if line matches pattern
    end # end of reading status file
    return @request_id
  end # end of request_id method

  def status
    return @status if(@status != "RUNNING" && @status != "INITIALIZING")
    status_file = "/tmp/deleteme.txt"
    command = "#{Dependencies[:html2text]} #{AServerConsole}/#{RequestString}#{@request_id} > #{status_file}"
    system(command)
    if($? == 0)
      pattern = /Current State: ([A-Z]+)/
      open(status_file).each do |line|
        if(line =~ pattern)
          @status = line.match(pattern).to_a[1]
          break
        end # end of check if line matches pattern
      end # end of reading status file

      # if there is was an assembly run-time error, copy log files
      if(@status == "ERROR" || @status == "CANCEL")
        error_dir = "#{@test_dir}/#{test_spec.version_label}_#{@genome.name}_error_logs"
        error_dir = Dir.create_path(error_dir)
        if(error_dir)
          command = "cp -Rp #{@aserver_dir}/command.out #{error_dir}"
          system(command)
          command = "cp -Rp #{@aserver_dir}/log #{error_dir}"
          system(command)
        end
      end
    else
      @status = "INDETERMINATE"
    end # end of if command ran ok
    @status
  end # end of status method

  def run_tampa
    return if(@status != "FINISHED")
    if(!File.exists?("#{@aserver_dir}/TAMPA"))
      Dir.chdir(@aserver_dir)
      command = "#{@test_spec.binary_path}/asm2TampaResults -a #{@genome.name} -b #{@test_spec.binary_path}"
      system(command)
    end

    # TAMPA has been run. Check that all output files are present
    @status = "TAMPA FINISHED"
    TampaFileSuffixes.values.each do |suffix|
      if(!File.exists?("#{@aserver_dir}/#{@genome.name}.#{suffix}"))
        @status = "TAMPA ERROR"
        break;
      end # check if file is present
    end # loop over TAMPA output file suffixes
  end # end of run_tampa method


end # end of assembler test class
######################################################################


######################################################################
# Genome Class
class Genome
  attr_reader :name, :dir, :ready

  def initialize(genome_name)
    @name = genome_name
    @dir = nil
    @ready = false
  end

  def set_up(base_dir, logger)
    return if(@ready)
    begin
      @dir = base_dir + "/" + @name
      @dir = Dir.create_path(@dir)
      if(!@dir)
        logger.log("Failed to create directory #{@dir}")
        return false
      end
      
      Dir.chdir(@dir)
      # check-for required genome files
      run_pullfrag = false
      InputFileSuffixes.each do |symb, suffix|
        fn = "#{@name}.#{suffix}"
        if(!File.exists?(fn))
          run_pullfrag = true
          break;
        end
      end
      
      # run pullfrag
      if(run_pullfrag)
        command = "#{Dependencies[:pullfrag]} -i -D #{@name} -o #{@name}"
        logger.log("Running #{Dependencies[:pullfrag]}", true)
        logger.log("\t#{command}")
        system(command)
        if($? != 0)
          # pullfrag failed
          logger.log("Failed to run #{Dependencies[:pullfrag]}", true)
          return false
        else
          # pullfrag ran to completion, but did it get what is needed?
          InputFileSuffixes.each do |symb, suffix|
            fn = "#{@name}.#{suffix}"
            if(!File.exists?(fn))
              logger.log("#{symb.to_s} file #{fn} not present", true)
              return false
            end
          end
        end
      end
      @ready = true
    ensure
      Dir.chdir(base_dir)
    end
    return true
  end
end # end of genome class
######################################################################


######################################################################
# "main" processing begins here

run_dir = File.original(Dir.pwd)

########################################
# process the command line
opts = GetoptLong.new(["--help", "-h", GetoptLong::NO_ARGUMENT],
                      ["--verbose", "-v", GetoptLong::NO_ARGUMENT],
                      ["--testspec", "-t", GetoptLong::REQUIRED_ARGUMENT],
                      ["--Genomes", "-G", GetoptLong::REQUIRED_ARGUMENT],
                      ["--genome", "-g", GetoptLong::REQUIRED_ARGUMENT],
                      ["--comp", "-c", GetoptLong::REQUIRED_ARGUMENT])

be_verbose = false
genome_files = {}
genomes = {}
test_spec_files = {}
comp_script_path = nil

opts.each do |opt, arg|
  case opt
  when "--help"
    print_help
  when "--verbose"
    be_verbose = true
  when "--testspec"
    fn = File.original(arg)
    print_help("Non-existent config file #{arg}") if(!fn)
    test_spec_files[fn] = (test_spec_files[fn] || test_spec_files.size + 1)
  when "--Genomes"
    fn = File.original(arg)
    print_help("Non-existent genomes file #{arg}") if(!fn)
    genome_files[fn] = true
  when "--genome"
    genome_name = arg.downcase
    genomes[genome_name] ||= Genome.new(genome_name)
  when "--comp"
    comp_script_path = File.original(arg)
  end
end

# override verbosity with testing
be_verbose = true if(TESTING || STUBBING)

# quick check of quickly checkable parameters
print_help("No test specification files") if(test_spec_files.size == 0)
print_help("No comparison script path") if(!comp_script_path)

# end of direct command line pricessing
########################################


########################################
# determine test number and log and results file names
# today is
StartDate = Time.now.strftime("%Y%m%d")

# determine test number needed for subdirectory naming
pattern = /^ca_test_#{StartDate}_([0-9]{3})$/
max_test_number = 0
Dir.foreach(run_dir) do |item_name|
  if(item_name =~ pattern)
    m = item_name.scan(pattern)[0][0].to_i
    max_test_number = (max_test_number < m ? m : max_test_number)
  end
end
TestNumber = sprintf("%03d", max_test_number + 1)
TestDir = "#{run_dir}/ca_test_#{StartDate}_#{TestNumber}"
if(!Dir.create_path(TestDir))
  printf(STDERR, "Failed to create test directiry #{TestDir}. Aborting.")
  exit 1
end

LogFile = "#{TestDir}/log.txt"
Log = Logger.new(LogFile, be_verbose)
Log.log("#{$0} started on #{ENV["HOSTNAME"]}.", true)

# end of log/result file name and test number determination
########################################


########################################
# load in genomes from genomes file(s)
genome_files.each_key do |fn|
  open(fn).each do |line|
    next if(line =~ /^#/)
    line.gsub!(/[\r\n\cZ]/, '')
    genome_name = genome.strip.split.compact[0].strip.downcase
    genomes[genome_name] ||= Genome.new(genome_name)
  end
end ## end of loop over genome files
if(genomes.size == 0)
  Log.log("No genomes specified.\nExiting.")
  exit 1
end
# end of genome file loading
########################################


########################################
# load in versions/run mods from config file(s)
# go through them in the order specified to determine baseline order

# keep track of what makes a test spec unique with a hash
test_specs = {}

# keep track of what test specs are labeled with a hash
version_labels = {}

# keep track of everything in baseline order
spec_order = {}

# iterate over config files
test_spec_files.sort{|a,b| a[1] <=> b[1]}.each do |fn, count|
  # iterate over lines in the config file
  open(fn).each do |spec_line|
    next if(spec_line =~ /^#/)

    # create a test spec object from the line
    test_spec = AssemblerTestSpec.new(spec_line)
    
    # check for spec errors
    errors = test_spec.errors
    if(errors)
      Log.log("\nErrors in test specification:")
      Log.log("  #{spec_line}")
      Log.log("#{errors}")
      next
    end

    # check for uniqueness of subdir names
    if(test_specs[test_spec.key] || version_labels[test_spec.version_label])
      Log.log("\nIgnoring duplicate test specification:")
      Log.log("#{test_spec.line}")
      Log.log("Version label:")
      Log.log("  #{test_spec.version_label}")
      Log.log("encountered previously in config file:")
      if(test_specs[test_spec.key])
        Log.log("  #{spec_order[test_specs[test_spec.key]][:config_file]}")
      else
        Log.log("  #{spec_order[version_labels[test_spec.version_label]][:config_file]}")
      end
      next
    else
      i = spec_order.size
      spec_order[i] = {}
      spec_order[i][:spec] = test_spec
      spec_order[i][:config_file] = fn
      test_specs[test_spec.key] = i
      version_labels[test_spec.version_label] = i
    end ## end of if item is in the test_specs hash
  end ## end of loop over lines in config file
end ## end of loop over config files

# check that one or more test specs were actually provided
if(spec_order.size == 0)
  Log.log("No test specifications.\nExiting.")
  exit 1
end
# end of config file loading
########################################


########################################
# Generate initial log entries reflecting inputs
Log.log("########################")
Log.log("# Command line feedback:")
Log.log("# Test directory:")
Log.log("#   #{TestDir}")
Log.log("# Genomes to assemble:")
Log.log("#   #{genomes.keys.sort.join(", ")}")
Log.log("# Test specifications in baseline order:")
spec_order.keys.sort.each do |i|
  Log.log("# #{i + 1}. #{spec_order[i][:spec].version_label}")
  Log.log("#       #{spec_order[i][:spec].key}")
end
Log.log("\n\n")
# end of initial log entries
########################################


########################################
# Run tests

# move to next phase - setting up runs
genomes.values.each do |genome|
  Dir.chdir(run_dir)
  
  if(!genome.set_up(run_dir, Log))
    Log.log("Failed to create genome #{genome.name} input files.")
    Log.log("Skipping #{genome.name}.")
    next
  end

  Dir.chdir(TestDir)
  assembler_tests = []
  
  # loop over test specs and execute them
  spec_order.keys.sort.each do |i|

    assembler_tests << AssemblerTest.new(spec_order[i][:spec],
                                         genome,
                                         TestDir)
    assembler_tests[-1].set_up
    if(assembler_tests[-1].status == "READY")
      assembler_tests[-1].run
      sleep(1)
      Log.log("#{spec_order[i][:spec].version_label} started on #{genome.name} with request id #{assembler_tests[-1].request_id}", true)
    end
  end

  # check/monitor runs
  completed_assemblies = {}
  while(completed_assemblies.size < assembler_tests.size)
    sleep(5)
    assembler_tests.each do |assembler_test|
      status = assembler_test.status
      case status
      when "INDETERMINATE"
        completed_assemblies[assembler_test] = true
      when "INITIATION ERROR"
        completed_assemblies[assembler_test] = true
      when "ERROR"
        completed_assemblies[assembler_test] = true
      when "CANCEL"
        completed_assemblies[assembler_test] = true
      when "FINISHED"
        completed_assemblies[assembler_test] = true
        assembler_test.run_tampa
      end
    end # loop over all tests to check status
  end # loop until all tests have finished

  # compare results
  # temporarily using qc_combine instead
  # qcCompare_command = "#{comp_script_path}/qcCompare"
  qcCompare_command = "qc_combine"
  tampaCompare_command = "#{comp_script_path}/tampaCompare"
  num_successful_qcs = 0
  num_successful_tampas = 0
  spec_order.keys.sort.each do |i|
    #qcCompare
    # rather than check for the status text, look for the .qc file
    qc_filename = assembler_tests[i].aserver_dir + "/" + genome.name + ".qc"
#    if(assembler_tests[i].status =~ /FINISHED/)
    if(File.exists?(qc_filename))
      num_successful_qcs += 1
      # qcCompare_command += " -f " + qc_filename
      qcCompare_command += " " + qc_filename
    else
      Log.log("Assembly #{assembler_tests[i].test_spec.version_label} failed for #{genome.name}")
    end

    #tampaCompare
    if(assembler_tests[i].status == "TAMPA FINISHED")
      num_successful_tampas += 1
      tampaCompare_command += " -d #{assembler_tests[i].aserver_dir} -a #{genome.name}"
    else
      Log.log("TAMPA failed to run on #{assembler_tests[i].test_spec.version_label} results for #{genome.name}")
    end
  end # loop over all tests for comparisons

  if(num_successful_qcs > 1)
    qcCompare_command += " > #{TestDir}/#{genome.name}_results.txt"
    system(qcCompare_command)
    if($? == 0)
      Log.log("Completed qcCompare on #{genome.name} tests", true)
    else
      Log.log("Failed to run qcCompare on #{genome.name} tests", true)
    end

    if(num_successful_tampas > 1)
      tampaCompare_command += " >> #{TestDir}/#{genome.name}_results.txt"
      system(tampaCompare_command)
      if($? == 0)
        Log.log("Completed tampaCompare on #{genome.name} tests", true)
      else
        Log.log("Failed to run tampaCompare on #{genome.name} tests", true)
      end
    else
      Log.log("TAMPA succeeded on fewer than 2 tests for #{genome.name}.")
      Log.log("No TAMPA comparisons to run.")
    end

  else
    Log.log("Fewer than 2 tests ran to completion for #{genome.name}.")
    Log.log("No QC or TAMPA comparisons to run.")
  end # check if there are two or more assemblies to compare
  
  Log.log("Done.", true)
end
