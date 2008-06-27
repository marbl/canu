#!/usr/bin/env ruby

asmFileS = ARGV[0]
ARGV.shift
asmFile = File.open( asmFileS )
mates = {}
ARGV.each do |frgFileS|
    frgFile = File.open( frgFileS )
    frgFile.each_line do |line|
        if line[0,4] == '{LKG'
            frgFile.readline
            frgFile.readline
            fg1L = frgFile.readline
            raise "Bad fg1 line: #{fg1L}" unless fg1L[0,4] == 'fg1:'
            fg2L = frgFile.readline
            raise "Bad fg2 line: #{fg2L}" unless fg2L[0,4] == 'fg2:'
            fg1L.chop!
            fg2L.chop!
            fg1 = fg1L[4,fg1L.length]
            fg2 = fg2L[4,fg2L.length]
            mates[ fg1 ] = fg2
            mates[ fg2 ] = fg1
        end
    end
end
t = '1096625850513'
puts "Mate of #{t} is #{mates[t]}" if mates.has_key?( t )

surrogates = {}
surrUTGs   = {}
placedSurrs= {}
allFrags   = {}
numAlreadySurr  = 0
numPlacedMate   = 0
numUnplacedMate = 0
numPlacedSurr   = 0
numUnplacedSurr = 0
asmFile.each_line do |line|
    if line[0,4] == '{AFG'
        line = asmFile.readline
        fid,rest = line[5,line.length].split(',')
        allFrags[ fid ] = 1

    elsif line[0,4] == '{UTG'
        acc,iid = asmFile.readline.chop[5..-1].split(',')
        line = asmFile.readline while line[0,4] != 'sta:'
        if line[4,1] == 'S'
            line = asmFile.readline while line[0,4] != 'nfr:'
            puts "Surrogate #{acc} #{line}"
            line.chop!
            surrUTGs[ acc ] = [ 0 ]
            nfr = line[4,line.length]
            nfr.to_i.times do
                line = asmFile.readline while line[0,4] != 'mid:'
                line.chop!
                mid = line[4,line.length]
                surrUTGs[ acc ].push( mid )
                if surrogates.has_key?( mid )
                    numAlreadySurr += 1
                else
                    surrogates[ mid ] = acc
#                        surrogates[ mates[ mid ] ] = 1
                end
                line = asmFile.readline
            end
        end
    elsif line[0,4] == '{CCO'
        line = asmFile.readline while line[0,4] != 'pla:'
        placed = line[4,1] == 'P'
        line = asmFile.readline while line[0,4] != 'npc:'
        nou  = asmFile.readline.chop[4..-1].to_i
#        puts "Placed? #{ placed }, num cco frags #{line}"
        line.chop!
        npc = line[4,line.length]
        npc.to_i.times do
            line = asmFile.readline while line[0,4] != 'mid:'
            line.chop!
            mid = line[4,line.length]
            if mates.has_key?( mid )
                if placed
                    numPlacedMate += 1
                else
                    numUnplacedMate += 1
                end
            end
            if surrogates.has_key?( mid )
                numAlreadySurr += 1
                if placed
                    numPlacedSurr += 1
                    placedSurrs[ mid ] = 1
                else
                    numUnplacedSurr += 1
                end
            end
            line = asmFile.readline
        end
        nou.times do
            line = asmFile.readline while line[0,4] != 'lid:'
            lid = line.chop[4..-1]
            if surrUTGs.has_key?( lid )
                surrUTGs[lid][0] += 1
            end
            line = asmFile.readline
        end
    end
end
allFrags.each_key do |frag|
    if mates.has_key?( frag ) && !allFrags.has_key?(mates[frag]) &&
       surrogates.has_key?( frag )
        puts "Surrogate frag #{frag} has mate, but it's not in .asm file."
    end
end
multiSurr = {}
surrUTGs.each_pair do |acc,frags|
    puts "Surrogate #{acc} occurs #{frags[0]} times."
    if frags[0] > 1
        frags.each { |frg| multiSurr[ frg ] = 1 unless frg == 0 }
    end
end
numBothSurr = numSurr = numMulti = numMultiBoth = numPlacedMulti = numPlacedBoth = 0
surrogates.each_key do |frag|
    if mates.has_key?( frag )
        if surrogates.has_key?( mates[frag] )
            numBothSurr += 1
            if multiSurr.has_key?( frag )
                numMultiBoth += 1
                numPlacedBoth += 1 if placedSurrs.has_key?( frag )
            end
        else
            numSurr += 1
            if multiSurr.has_key?( frag )
                numMulti+= 1
                numPlacedMulti += 1 if placedSurrs.has_key?( frag )
            end
        end
    end
end
puts "Num surrogate reads #{ surrogates.length }"
puts "Num surrogate mates #{ numSurr }"
puts "Num both surrogate mates #{ numBothSurr }"
puts "Num multi surrogate mates #{ numMulti }"
puts "Num multi both surrogate mates #{ numMultiBoth }"
puts "Num multi surrogate mates placed #{ numPlacedMulti } unplaced #{numMulti - numPlacedMulti}"
puts "Num multi both surrogate mates placed #{ numPlacedBoth } unplaced #{numMultiBoth - numPlacedBoth}"
puts "Num surrogate already appeard #{ numAlreadySurr }"
puts "Num placed mated reads #{ numPlacedMate }"
puts "Num unplaced mated reads #{ numUnplacedMate }"
puts "Num placed mated surrogate reads #{ numPlacedSurr }"
puts "Num unplaced mated surrogate reads #{ numUnplacedSurr }"
