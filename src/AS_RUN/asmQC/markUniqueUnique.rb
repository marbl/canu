#!/usr/bin/env ruby

asmFileS = ARGV[0]
iumFileS = ARGV[1]

asmFile = File.open( asmFileS )
iumFile = File.open( iumFileS )


surrUTGs   = {}
iidToUid   = {}
placedSurrs= {}
asmFile.each_line do |line|

    if line[0,4] == '{UTG'
        acc,iid = asmFile.readline.chop[5..-1].split(',')
        line = asmFile.readline while line[0,4] != 'sta:'
        if line[4,1] == 'S'
            iidToUid[ iid.chop ] = acc
            line = asmFile.readline while line[0,4] != 'nfr:'
            #            puts "Surrogate #{acc} #{line}"
            line.chop!
            surrUTGs[ acc ] = [ 0 ]
            nfr = line[4,line.length]
            nfr.to_i.times do
                line = asmFile.readline while line[0,4] != 'mid:'
                line.chop!
                mid = line[4,line.length]
                surrUTGs[ acc ].push( mid )
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
uniqueSurr = {}
surrUTGs.each_pair do |acc,frags|
    $stderr.puts "Surrogate #{acc} occurs #{frags[0]} times."
    if frags[0] == 1
        uniqueSurr[ acc ] = 1
    end
end

iidToUid.each_pair { |iid,uid| $stderr.puts "iid #{iid} uid #{uid}" }

# read and output IUM file changing unique surrogate status to unique
iumFile.each_line do |line|

    if line[0,4] == '{IUM'
        print line
        line = iumFile.readline
        accStr,acc = line.chop.split(':')
        uid = iidToUid[ acc ]
        if uniqueSurr[ uid ] == 1
            # read up to fur: line
            while line[0,4] != 'fur:' do
                print line
                line = iumFile.readline
            end
        
            uniqueStat = line[4,1]
            if uniqueStat != 'U'
                line = "fur:U\n"
                $stderr.print "Marking unitig #{uid},#{acc} unique was #{uniqueStat}\n"
            end
        end
    end
    # echo back most lines
    print line
end
