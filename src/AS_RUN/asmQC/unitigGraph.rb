#!/usr/bin/env ruby

require 'rgl/adjacency'
require 'rgl/traversal'
require 'rgl/implicit'
require 'rgl/mutable'
require 'rgl/rdot'
#require 'graph/graphviz_dot'

usage = "#$0 unitig_iid asm/iumFile [UTG|CCO] [bestEdgeFile]
Give it an iid for either a contig(CCO) or a unitig(UTG),
a file with messages to parse and optionally
tell it if you want unitigs or contigs. Default is contigs.

Currently needs env var RUBYLIB=~eventer/lib/ruby for rgl module.
"

node = ARGV[0]
iumFilePath = ARGV[1]
utgStr = ARGV[2]
bestEdgeFilePath = ARGV[3]
abort usage unless node && iumFilePath

# default to contig message, just because it's smaller
$tigMsg  = '{CCO'
$linkMsg = '{CLK'
if utgStr == 'UTG'
    $tigMsg  = '{UTG'
    $linkMsg = '{ULK'
elsif utgStr && utgStr != 'CCO'
    abort usage
end

#iumFile      = File.open( iumFilePath )
asmFile      = File.open( iumFilePath )
unitigs = {}
#bestEdgeFile = File.open( bestEdgeFilePath )

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

class File
    def readTo(str)
        line = ''
        line = readline while line[0,4] != str
        line
    end
end

$edgeInfo = {}
$firstLast = {}
$fragsUnitig = {}
def readUnitigsFromIUMFile(iumFile)
    unitigs = {}
    iumFile.each_line do |line|
        if line[0,4] == '{IUM'
            ac,acc = iumFile.readline.chomp.split(':')
            raise "Expected acc:<id>" unless ac == 'acc'

            line = iumFile.readTo('nfr:')
            nfr = line.chop[4,line.length].to_i
            next if nfr == 1

            line = iumFile.readTo('mid:')
            mid = line.chop[4,line.length]
            $fragsUnitig[ mid ] = acc

#            next if nfr == 1
            utg = Unitig.new(acc)
            utg.firstFrag = mid

            line = iumFile.readTo('pos:')
            b,e = line.chop[4,line.length].split(',')
            whichEnd = "5'"
            whichEnd = "3'" if e < b

            $firstLast[ mid ] = whichEnd

            nfr -= 1
            lastNonContain = 0
            nfr.times do
                line = iumFile.readTo('mid:')
                mid = line.chop[4,line.length]
                utg.addInternalFrag( mid )
                $fragsUnitig[ mid ] = acc
                line = iumFile.readTo('con:')
                if line == "con:0\n"
                    lastNonContain = mid
                    b,e = iumFile.readline.chop[4..-1].split(',')
                    whichEnd = "3'"
                    whichEnd = "5'" if e < b
                end
            end
            if lastNonContain == 0
                $fragsUnitig.delete( utg.firstFrag )
                next
            end
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
$iidToLen = {}
$surrogates = {}
def readUnitigsFromAsmFile(asmFile)
    uidToIID = {}
    links = []
    ulk = {}
    asmFile.each_line do |line|
        if line[0,4] == $tigMsg
            uid,iid = asmFile.readline.scan(/\d+/)
            raise "Bad acc: #{uid}" if uid == nil || iid == nil
            if $tigMsg == '{UTG'
                surr = asmFile.readTo('sta:')[4]
                if surr.chr == 'S'
                    $surrogates[ iid ] = true
                end
            end
            line = asmFile.readTo('len:')
            len = line[4,line.length].chomp
#            next unless len.to_i > 1000
            uidToIID[ uid ] = iid
            $iidToLen[ iid ] = len

        elsif line[0,4] == $linkMsg
            u,ut1 = asmFile.readline.scan(/\d+/)
            u,ut2 = asmFile.readline.scan(/\d+/)
            ut1 = uidToIID[ ut1 ]
            ut2 = uidToIID[ ut2 ]
            next unless ut1 != nil && ut2 != nil
            line = asmFile.readline
            ori = line[4]
            o1= o2 = color =''
            case ori
            when 78 then o1,o2,color = "3'","5'",'blue'
            when 65 then o1,o2,color = "5'","3'",'green'
            when 79 then o1,o2,color = "5'","5'",'yellow'
            when 73 then o1,o2,color = "3'","3'",'red'
            else raise "Invalid ori: #{ori}, #{line}"
            end
            line = asmFile.readTo('mea:')
            mean = line[4,line.length].to_i
            line = asmFile.readTo('num:')
            num = line[4,line.length].to_i
#            next if num < 2 || mean < 1
#            next if num < 2
            if ulk.has_key?( ut1 )
                if ulk[ut1].has_key?(ut2)
                    ulk[ut1][ut2] += 1
                else
                    ulk[ut1][ut2] = 1
#                    link = Link.new(ut1, ut2, "#{ut1} #{ut2}")
                    $edgeInfo["#{ut1} #{ut2}"] = "#{color},#{o1}  #{mean}  #{o2}"
                    links.push( [ut1,ut2] )
                end
            else
                ulk[ut1] = { ut2 => 1 }
#                link = Link.new(ut1, ut2, "#{ut1} #{ut2}")
                $edgeInfo["#{ut1} #{ut2}"] = "#{color},#{o1}  #{mean}  #{o2}"
                links.push( [ut1,ut2] )
            end
        end
    end
    return links,ulk
end

def graphrStuff()
    dgp = DotGraphPrinter.new(links)
    dgp.node_shaper = proc {|n| 'circle'}
#dgp.link_labeler = proc { |info| ["\"#{info[0,2]}\"","\"#{info[3,info.length]}\""] }
    dgp.link_labeler = proc { |info|
        u1,u2 = info.split
        cnt = ulk[u1][u2]
        color = 'black'
        if 15 < cnt
            color = 'yellow'
        elsif 11 < cnt
            color = 'orange'
        elsif 7 < cnt
            color = 'red'
        elsif 4 < cnt
            color = 'green'
        elsif 2 < cnt
            color = 'blue'
        end
        ["\"#{ulk[u1][u2]}\"","\"#{color}\""]
    }
#puts dgp.to_dot_specification
    dgp.write_to_file("utggraph.png","png")
    dgp.orientation = "landscape"      # Dot problem with PS orientation
    dgp.write_to_file("utggraph.ps")          # Generate postscript file
end

#unitigs = readUnitigsFromIUMFile( iumFile )

#bestEdges = readBestEdgeFile( bestEdgeFile )
#puts "Best edge for 15824 is #{bestEdges['15824']}, 13920 is #{bestEdges['13920']}"
#puts "Best edge for 27 is #{bestEdges['27']}, 36838 is #{bestEdges['36838']}"
#puts "Best edge for 12124 is #{bestEdges['12124']}, 36480 is #{bestEdges['36480']}"
graph = RGL::AdjacencyGraph.new()
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
    elsif not $fragsUnitig.has_key?(lastEdge)
        puts "Skip #{acc} last frag #{utg.lastFrag} singleton edge."
    else
#        graph.add_edge( utg.accession, $fragsUnitig[lastEdge])
        color = 'black'
        otherUtg = $fragsUnitig[lastEdge]
        raise "Edge #{lastEdge} from #{acc} to bad untig #{otherUtg}" if not
         unitigs.has_key?( otherUtg )
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
def findAdjacentAtDepth(graph, startNode, maxDepth)
    links = []
    used = {}
    nodes = [startNode]
    adj = []
    (maxDepth+1).times {
        nodes.each { |node|
            if used.has_key?( node )
                next
            else
                used[ node  ] = true
            end
            adj.concat( graph.adjacent_vertices( node ) )
        }
        nodes = Array.new(adj)
        adj.clear
    }
    graph.each_edge { |u,v|
        if used.has_key?( u ) && used.has_key?( v )
                links.push( [u, v] )
        end
    }
    links
end
def graphToRDot(graph)
    nodes = {}
    graph.each_edge { |u,v|
        key = "#{u} #{v}"
        yek = "#{v} #{u}"
        key,u,v = yek,v,u if $edgeInfo.has_key?( yek )
        if $edgeInfo.has_key?( key )
            [u,v].each { |l| lab = "\"#{l} :#{$iidToLen[l]}\""
                unless nodes.has_key?(l)
                    shape = 'ellipse'
                    shape = 'triangle' if $surrogates[ l ]
                    nodes[l]=DOT::DOTNode.new({'name' => l,'label' => lab,
                                        'fontsize' => 8, 'shape' => shape
                    #                    , 'height' => 0.2, 'width' => 0.5
                    })
                end
            }
            color,label = $edgeInfo[key].split(',')

            edge = DOT::DOTDirectedEdge.new( {"from" => u, "to" => v,
                        "label" => label,
                        "color" => color,
                        "fontsize" => 7}
            )
            nodes[ u ] << edge
        end
    }
    DOT::DOTDigraph.new({'name' => 'UnitigGraph','nodes' => nodes.values})
end
class DOT::DOTSubgraph
    def write_rdot_to_graphic_file (fmt='png', dotfile="graph")
        src = dotfile + ".dot"
        dot = dotfile + "." + fmt

        File.open(src, 'w') do |f|
            f << self.to_s << "\n"
        end

        system( "dot -T#{fmt} #{src} -o #{dot}" )
        dot
    end
end

links,ulk = readUnitigsFromAsmFile(asmFile)
graph.add_edges(*links)
#visitor = RGL::DFSVisitor.new(graph)
#visitor.attach_distance_map
#node = '235'
#node = '67429'
#node = '1097326128195'
#node = '0'
#graph.write_to_graphic_file('ps','allutg')
#graph.each_vertex {|v| puts "#{v}"}
#puts graph.adjacent_vertices( node )
#graph.depth_first_visit( node, visitor) {}
#fgraph = graph.vertices_filtered_by {|v| visitor.distance_to_root(v) < 4}
#fgraph = findAdjacentAtDepth(graph, node, 2)
links = findAdjacentAtDepth(graph, node, 20)
#fgraph.each_edge { |u,v| puts "#{u}-#{v}" }
fgraph = RGL::AdjacencyGraph.new()
fgraph.add_edges(*links)
#fgraph.write_to_graphic_file('ps',"unitigOrig#{node}")
dotGraph = graphToRDot( fgraph )
#puts dotGraph
dotGraph.write_rdot_to_graphic_file('ps',"#{utgStr}#{node}")
