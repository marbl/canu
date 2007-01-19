#!/usr/bin/env ruby

gkpStore = (ARGV[0] || '*.gkpStore')
bin = ''
if ARGV[1]
    bin = ARGV[1]
    bin +='/' unless bin[-1] == '/'
end

cmd = %Q[#{bin}dumpGatekeeper #{gkpStore} 2>/dev/null | grep 'Link (']
dump = IO.popen(cmd)
raise "popen failure: #{cmd}" if dump.eof

re = /^\tLink \((\d+),(\d+)\) dist:(\d+) /
ranges = {}
dump.each_line do |line|
   frg1,frg2,dist = re.match(line)[1,3].collect{|s| s.to_i}
   if frg1 > frg2 then raise "Order mismatch #{frg1} > #{frg2}" end
   if ranges.has_key?( dist )
     if ranges[dist][1] < frg2
        ranges[dist][1] = frg2
     end
   else
     ranges[dist] = [frg1,frg2]
   end
end
ranges.each {|key,val| puts "Distance #{key} has IID range #{val[0]} to #{val[1]}"}
ranges.each {|key1,val1|
    ranges.each {|key2,val2|
        next if key1 == key2
        if val1[0].between?(val2[0],val2[1]) ||
           val1[1].between?(val2[0],val2[1]) 
           raise "Range overlap error Distance #{key1} and #{key2} overlap"
        end
    }
}
