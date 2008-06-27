#!/usr/bin/ruby

class Scaffold
    attr_accessor :bases, :gapsize, :contigs
    def initialize()
        @bases = 0
        @gapsize = 0
        @contigs = 0
    end
end

scafs = {}
totalBases = 0
totalGaps = 0
biggest = 0
$stdin.each_line do |line|
    scf,cbeg,cend,cnum,typ,ctg,one,clen,ort = line.chomp.split
    if typ == 'W'
        if scafs.has_key?(scf)
            sObj = scafs[scf]
        else
            sObj = Scaffold.new
            scafs[scf] = sObj
        end
        sObj.bases += clen.to_i
        totalBases += clen.to_i
        sObj.contigs += 1
        biggest = sObj.bases if sObj.bases > biggest
    elsif typ == 'N'
        totalGaps += ctg.to_i
        scafs[ scf ].gapsize += ctg.to_i
    else
        raise "Unknown Line in AGP file: #{ line }"
    end
end

puts "Total bases in scaffolds is #{totalBases}"
puts "Total bases in scaffold gaps is #{totalGaps}"
puts "The longest scaffold has #{biggest} bases"
numScafs = scafs.length
puts "There are #{numScafs} scaffolds, for an avg size of #{totalBases/numScafs}"
