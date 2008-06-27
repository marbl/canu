#!/usr/bin/env ruby

class File
    def readTo(str)
        line = ''
        line = readline while line[0,4] != str
        line
    end
end

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
        line = asmFile.readTo( 'sta:' )
        if line[4,1] == 'S' || line[4,1] == 'N'
            iidToUid[ iid.chop ] = acc
            surrLen = asmFile.readTo( 'len:' ).chop[4..-1]
            next if surrLen.to_i < 2000
            nfr = asmFile.readTo( 'nfr:' ).chop[4..-1]
            surrUTGs[ acc ] = [ 0 ]
            nfr.to_i.times do
                mid = asmFile.readTo( 'mid:' ).chop[4..-1]
                surrUTGs[ acc ].push( mid )
            end
        end
    elsif line[0,4] == '{CCO'
        line = asmFile.readTo( 'pla:' )
        placed = line[4,1] == 'P'
        npc = asmFile.readTo( 'npc:' ).chop[4..-1]
        nou  = asmFile.readline.chop[4..-1].to_i
#        puts "Placed? #{ placed }, num cco frags #{line}"
        npc.to_i.times do
            mid = asmFile.readTo( 'mid:' )
        end
        nou.times do
            lid = asmFile.readTo( 'lid:' ).chop[4..-1]
            if surrUTGs.has_key?( lid )
                surrUTGs[lid][0] += 1
            end
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

#iidToUid.each_pair { |iid,uid| $stderr.puts "iid #{iid} uid #{uid}" }

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
