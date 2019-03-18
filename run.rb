#!/usr/bin/ruby
require 'optparse'
require 'fileutils'
require '/home/cew54/slurmer/Qsub.rb'
require 'pp'

load("~/DIRS.txt") # gets MIKHAILDIR

OPTIONS = {}
OPTIONS[:int] = false
OPTIONS[:run] = false
OPTIONS[:norun] = false
OptionParser.new do |opts|
  opts.banner = "Usage: run.rb [OPTIONS] COMMAND"

  opts.on("-i", "--[no-]interactive", "Run each job in serial on the interactive log in node") do |i|
    OPTIONS[:int] = i
  end
  
  opts.on("-n", "--norun", "Don't run any jobs, just print running messages") do |i|
    OPTIONS[:norun] = i
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
           coloc :: run coloc for a subset of things

           clean :: remove slurm files      
  "
end
if OPTIONS[:help] then
  usage()
  exit 0
end


comfile="run-#{COMMAND}.sh"
args = {:job=>COMMAND,
        :tasks=>'1',
        :cpus=>'1',
:time=>'04:00:00',
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
  puts "--- prep-files-v5.R ---"
  commands[0] = "PROX=0 Rscript ./prep-files-v5.R > log/prep-files-v5-0.log"
  commands[1] = "PROX=1 Rscript ./prep-files-v5.R > log/prep-files-v5-1.log"
  # commands[0] = "Rscript ./prep-files-v4.R --args PROP=TRUE > log/prep-files-v4.log"

when "guess"
  puts "--- GUESS ---"
  args = {:job=>COMMAND,
        :tasks=>'1',
        :cpus=>'1',
        :autorun=>OPTIONS[:autorun],
        :time=>'08:00:00',
        :excl=>" "}
  comfiles=listfiles(MIKHAILDIR, "results/*/runme.sh")
  resfiles=listfiles(MIKHAILDIR, "results/*/out_500000_features.txt")
  skipfiles=listfiles(MIKHAILDIR, "results/*/skip")
  todo=comfiles.map { |f| File.dirname(f) } - resfiles.map { |f| File.dirname(f) } - skipfiles.map { |f| File.dirname(f) }
  todo=todo.map { |f| f + "/runme.sh" }
  puts "guess, wanted: #{comfiles.length}"
  puts "guess, done: #{resfiles.length}"
  puts "guess, skip: #{skipfiles.length}"
  puts "guess, todo: #{todo.length}"
  commands=todo.map { |f| File.read(f).chomp }

when "expand"
  puts "--- EXPAND ---"
  ffiles=listfiles(MIKHAILDIR,"results/*/out_500000_features.txt")
  done=listfiles(MIKHAILDIR,"results/*/*-nsnp.csv")
  todo = ffiles.map { |f| File.dirname(f) } - done.map { |f| File.dirname(f) }
  puts "expand, done: #{done.length}"
  puts "expand, todo: #{todo.length}"
  commands = todo.map { |f| "./expand-results.R --args d=" + f }

when "coloc"
  puts "--- COLOC ---"
    args = {:job=>COMMAND,
        :tasks=>'1',
        :cpus=>'1',
        :autorun=>OPTIONS[:autorun],
        :time=>'1:00:00',
        :excl=>" "}
    # com="awk '$3==\"TRUE\" && $14==\"TRUE\" {print $1,$15}' #{MIKHAILDIR}/all_tested_affinity_expression_associations_v3.txt"
    # com="awk '$3==\"TRUE\" && $6==\"TRUE\" {print $1,$7}' #{MIKHAILDIR}/all_tested_affinity_expression_associations_v5.txt"
    com="awk -F, '{print $2,$3}' assoc-prox.tab"
    outfiles = listfiles(MIKHAILDIR, "coloc-v5/*.csv")
    infiles = listfiles(MIKHAILDIR, "results/*/*-nsnp.csv").map{ |f| File.basename(File.dirname(f)) }

    genes=%x[#{com}].split("\n")
    genes.each { |gp|
      g=gp.split(" ")
      next unless infiles.include? g[0] and infiles.include? g[1]
      of="#{MIKHAILDIR}/coloc-v5/#{g[0]}-#{g[1]}.csv"
      commands.push("./coloc-v5.R --args g1=#{g[0]} g2=#{g[1]}") unless outfiles.include? of
    }
    puts "coloc todo: #{commands.length}"
# ENSG00000003402 ENSG00000240344
# ENSG00000003402 ENSG00000115942
# ENSG00000003402 ENSG00000183308
# ENSG00000003402 ENSG00000064012

when "colocgz"
   system("cd #{MIKHAILDIR} && cat coloc-v5/*.csv | awk 'NR==1 || $1!=\"nsnps\"' | gzip -c > coloc.gz")

when "tgz"
   system("cd #{MIKHAILDIR} && tar zcvf results.tgz results/*/*snpmod-99.RData results/*/*-nsnp.csv results/*/*-mppi.csv */*/skip */*/qcflag */*/ess.png")
end
################################################################################

## run

if(commands.length > 0 and !OPTIONS[:norun]) then
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




