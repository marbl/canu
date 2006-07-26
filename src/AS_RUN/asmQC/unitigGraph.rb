#!/usr/bin/env ruby

#require 'rgl/adjacency'
#require 'rgl/dot'
require 'graph/graphviz_dot'

class Link
    attr_accessor :from, :to, :info
    def initialize(from, to, info)
        @from = from
        @to = to
        @info = info
    end
end
class Unitig
    attr_accessor :firstFrag, :lastFrag
    attr_reader   :accession

    def initialize(acc)
        @accession = acc
        @firstFrag = 0
        @lastFrag = 0
        @internalFrags = {}
    end
    
    def has_frag?(frag)
        if @firstFrag == frag || @lastFrag == frag ||
           @internalFrags.has_key?(frag)
            return true
        else
            return false
        end
    end

    def addInternalFrag(frag)
        if has_frag?( frag )
            raise "Utg #{self.accession} already has frag #{frag}"
        end
        @internalFrags[ frag ] = 1
    end

    def numInternalFrags()
        @internalFrags.length
    end
    def deleteInternalFrag(frag)
        @internalFrags.delete(frag)
    end
end

$firstLast = {}
$fragsUnitig = {}
def readUnitigsFromIUMFile(iumFile)
    unitigs = {}
    iumFile.each_line do |line|
        if line[0,4] == '{IUM'
            ac,acc = iumFile.readline.chomp.split(':')
            raise "Expected acc:<id>" unless ac == 'acc'

            line = iumFile.readline while line[0,4] != 'nfr:'
            nfr = line.chop[4,line.length].to_i
            next if nfr == 1

            line = iumFile.readline while line[0,4] != 'mid:'
            mid = line.chop[4,line.length]
            $fragsUnitig[ mid ] = acc

#            next if nfr == 1
            utg = Unitig.new(acc)
            utg.firstFrag = mid
            
            line = iumFile.readline while line[0,4] != 'pos:'
            b,e = line.chop[4,line.length].split(',')
            whichEnd = "5'"
            whichEnd = "3'" if e < b 

            $firstLast[ mid ] = whichEnd

            nfr -= 1
            lastNonContain = 0
            nfr.times do
                line = iumFile.readline while line[0,4] != 'mid:'
                mid = line.chop[4,line.length]
                utg.addInternalFrag( mid )
                $fragsUnitig[ mid ] = acc
                line = iumFile.readline while line[0,4] != 'con:'
                if line == "con:0\n"
                    lastNonContain = mid
                    b,e = iumFile.readline.chop[4..-1].split(',')
                    whichEnd = "3'"
                    whichEnd = "5'" if e < b 
                end
            end
            next if lastNonContain == 0
            utg.deleteInternalFrag( lastNonContain )
            utg.lastFrag = lastNonContain
            $firstLast[ lastNonContain ] = whichEnd
            unitigs[acc] = utg
        end
    end
    return unitigs
end

def readBestEdgeFile( bestEdgeFile )
    bestEdge = {}
    bestEdgeFile.each_line do |line|
        frg,fv,c1,c2,c3,st1,f1,f2,f3,th,tc1,tc2,tc3,st2,t1,t2,t3 = line.chop.split
        next unless $firstLast.has_key?( frg )

        if st1 == 'best'
            th,tc1,tc2,tc3,st2,t1,t2,t3 = f2,f3,th,tc1,tc2,tc3,st2,t1
            f2 = f1
            f1 = f3 = 0
        end
        if st2 == 'best'
            t2 = t1
            t2 = t3 = 0
        end
        if $firstLast[ frg ] == fv
            bestEdge[frg] = f2
        elsif $firstLast[ frg ] == th
            bestEdge[frg] = t2
        else
            raise "Bad edge: #{frg}, #{$firstLast[frg]}"
        end
    end
    return bestEdge
end
iumFilePath      = ARGV[0]
bestEdgeFilePath = ARGV[1]
iumFile      = File.open( iumFilePath )
bestEdgeFile = File.open( bestEdgeFilePath )

unitigs = readUnitigsFromIUMFile( iumFile )

bestEdges = readBestEdgeFile( bestEdgeFile )
puts "Best edge for 15824 is #{bestEdges['15824']}, 13920 is #{bestEdges['13920']}"
puts "Best edge for 27 is #{bestEdges['27']}, 36838 is #{bestEdges['36838']}"
puts "Best edge for 12124 is #{bestEdges['12124']}, 36480 is #{bestEdges['36480']}"
#graph = RGL::DirectedAdjacencyGraph.new()
links = []
unitigs.each_pair do |acc,utg|
    firstEdge = bestEdges[ utg.firstFrag ]
    if firstEdge == nil
        raise "nil firstEdge"
    elsif firstEdge == '0'
        puts "Skip #{acc} first frag #{utg.firstFrag} with '0' edge."
    elsif firstEdge == 0
        puts "Skip #{acc} first frag #{utg.firstFrag} with 0 edge."
    elsif !$fragsUnitig.has_key?(firstEdge)
        puts "Skip #{acc} first frag #{utg.firstFrag} singleton edge."
    else 
#        graph.add_edge( utg.accession, $fragsUnitig[firstEdge])
        color = 'black'
        otherUtg = $fragsUnitig[firstEdge]
        if unitigs[otherUtg].firstFrag == firstEdge
            # begin to begin, make it green
            color = 'green'
        elsif unitigs[otherUtg].lastFrag == firstEdge
            # begin to end, make it red
            color = 'red'
        end
        link = Link.new(acc, otherUtg, "5' #{color}")
#links.push([acc, otherUtg, "5' color=#{color}"])
        links.push(link)
    end

    lastEdge  = bestEdges[ utg.lastFrag ]
    if lastEdge == nil
        raise "nil lastEdge #{acc} #{utg.lastFrag}"
    elsif lastEdge == '0'
        puts "Skip #{acc} last frag #{utg.lastFrag} with '0' edge."
    elsif lastEdge == 0 
        puts "Skip #{acc} last frag #{utg.lastFrag} with 0 edge."
    elsif !$fragsUnitig.has_key?(lastEdge)
        puts "Skip #{acc} last frag #{utg.lastFrag} singleton edge."
    else 
#        graph.add_edge( utg.accession, $fragsUnitig[lastEdge])
        color = 'black'
        otherUtg = $fragsUnitig[lastEdge]
        if unitigs[otherUtg].firstFrag == lastEdge
            # end to begin, make it green
            color = 'green'
        elsif unitigs[otherUtg].lastFrag == lastEdge
            # end to end, make it red
            color = 'red'
        end
        link = Link.new(acc, otherUtg, "3' #{color}")
#        links.push([acc, otherUtg, "3' #{color}"])
        links.push(link)
    end
end
dgp = DotGraphPrinter.new(links)
dgp.node_shaper = proc {|n| 'circle'}
dgp.link_labeler = proc { |info| ["\"#{info[0,2]}\"","\"#{info[3,info.length]}\""] }
#puts dgp.to_dot_specification
dgp.write_to_file("utggraph.png","png")
dgp.orientation = "landscape"      # Dot problem with PS orientation
dgp.write_to_file("utggraph.ps")          # Generate postscript file

#graph.write_to_graphic_file
