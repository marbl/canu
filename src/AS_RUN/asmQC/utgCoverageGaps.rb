#!/usr/bin/ruby

uidListFile = File.open(ARGV[0])
asmFile     = File.open(ARGV[1])
frgs = {}
uidListFile.each_line { |line| frgs[ line.chop ] = 1 }
asmFile.each_line do |line|
    if line[0,4] == '{UTG'
        line = asmFile.readline while line[0,4] != 'acc:'
        acc = line[4,line.length].chop

        line = asmFile.readline while line[0,4] != 'len:'
        len = line[4,line.length].to_i

        line = asmFile.readline while line[0,4] != 'cns:'
        cns = ''
        while line[0,1] != '.'
            line = asmFile.readline
            cns += line.chop
        end

        line = asmFile.readline while line[0,4] != 'nfr:'
        nfr = line[4,line.length].to_i
        mps = 0
        range = []
        while mps < nfr
            mps += 1
            line = asmFile.readline while line[0,4] != 'mid:'
            mid = line[4,line.length].chop
#puts "Looking at mid #{mid}"
            line = asmFile.readline
            next unless frgs.has_key?( mid )
#            puts "In frg list #{mid}"
            line = asmFile.readline while line[0,4] != 'pos:'
            b,e = line[4,line.length].chop.split(',')
            b,e = b.to_i,e.to_i
            b,e = e,b if b > e
            last = range.size - 1
            if last == -1 || b > range[last][1] # new range
                range.push([b,e])
            elsif e > range[last][1]
                range[last][1] = e
            end
        end
        $stderr.puts "UTG #{acc} of length #{len} is covered by frags in list from:"
        lastEnd = 0
        range.each_with_index do |r,idx|
           $stderr.puts "#{r[0]}-#{r[1]}"
           if r[0] - lastEnd > 80
               puts ">#{acc.delete('()').tr(',','|')}|#{idx} /offset=#{lastEnd}-#{r[0]}"
               puts cns[ lastEnd, r[0]-lastEnd ].delete('-')
           end
           lastEnd = r[1]
        end
    end
end
