#!/usr/bin/ruby
require 'optparse'
require 'fileutils'
require '/home/cew54/slurmer/Qsub.rb'
require 'pp'

load("~/DIRS.txt") # gets MIKHAILDIR

OPTIONS = {}
OPTIONS[:int] = false
OPTIONS[:run] = false
OptionParser.new do |opts|
  opts.banner = "Usage: runguess.rb [OPTIONS] COMMAND"

  opts.on("-i", "--[no-]interactive", "Run each job in serial on the interactive log in node") do |i|
    OPTIONS[:int] = i
  end
  
  opts.on("-r", "--autoRun", "Run each job on the queue without user input") do |i|
    OPTIONS[:autorun] = i
  end
  
  opts.on("-v", "--[no-]verbose", "Run verbosely") do |v|
    OPTIONS[:verbose] = v
  end

  opts.on("-h", "--help", "Show help") do |h|
    OPTIONS[:help] = h
  end
  
end.parse!
COMMAND = ARGV.shift

def usage()
  # puts OPTIONS
  puts "Usage: run.rb [OPTIONS] COMMAND

  COMMANDS are:
      (run things)
           prep :: prepare GUESS input
           guess :: run guess
           expand :: expand guess output

           clean :: remove slurm files      
  "
end
if OPTIONS[:help] then
  usage()
  exit 0
end


comfile="runguess-#{COMMAND}.sh"
args = {:job=>COMMAND,
        :tasks=>'1',
        :cpus=>'1',
           :autorun=>OPTIONS[:autorun],
           :excl=>" "}

def listfiles(dir,patt)
  files = Dir.glob(dir + '/' + patt)
  puts "FILES in #{dir}" if OPTIONS[:verbose]
  puts files if OPTIONS[:verbose]
  return files
end

def command_R(region,com,args='',parent='')
  if(parent=='..') then
    "./#{com}.R --args d=#{MIKHAILDIR}/#{region} #{args} > log/#{com}-#{region}.Rout 2>&1"
  else 
    "./#{com}.R --args d=#{MIKHAILDIR}/#{region}/GUESS #{args} > log/#{com}-#{region}.Rout 2>&1"
  end
end

## command
commands = []

case COMMAND
when "clean"
  rmfiles=listfiles(".", "runguess-*.sh*") +
          listfiles(".", "slurm-*.sh*") +
          listfiles(".", "machine.file.*")
  if(rmfiles.length() > 0) then
    puts "removing #{rmfiles.length()} files"
    puts(rmfiles)
    rmfiles.map { |f| File::delete(f) }
  end
  
when "prep"
  puts "--- prep-files.R ---"
  commands[0] = "Rscript ./prep-files.R > log/prep-files.log"

when "guess"
  puts "--- GUESS ---"
  args = {:job=>COMMAND,
        :tasks=>'1',
        :cpus=>'1',
        :autorun=>OPTIONS[:autorun],
        :time=>'04:00:00',
        :excl=>" "}
  comfiles=listfiles(MIKHAILDIR, "results/*/runme.sh")
  resfiles=listfiles(MIKHAILDIR, "results/*/out_500000_features.txt")
  todo=comfiles.map { |f| File.dirname(f) } - resfiles.map { |f| File.dirname(f) }
  todo=todo.map { |f| f + "/runme.sh" }
  commands=todo.map { |f| File.read(f).chomp }

when "expand"
  puts "--- EXPAND ---"
  ffiles=listfiles(MIKHAILDIR,"results/*/out_500000_features.txt")
  done=listfiles(MIKHAILDIR,"results/*/*-nsnp.csv")
  todo = ffiles.map { |f| File.dirname(f) } - done.map { |f| File.dirname(f) }
  commands = todo.map { |f| "./expand-results.R --args d=" + f }
end

################################################################################

## run

if(commands.length > 0) then
  if OPTIONS[:int] then
    puts "RUNNING INTERACTIVELY"
    commands.each { |s|
      puts s
      system(s)
    }
  else
    nodes = (commands.length / (args[:tasks].to_f)).ceil
    ncomm = commands.length
    puts "  RUNNING #{ncomm} commands over #{nodes} nodes"
    q=Qsub.new("runguess-#{COMMAND}.sh", args)
    commands.each { |s| q.add(s) }
    q.close()
  end
else
  puts "  NO COMMANDS TO RUN"
end




