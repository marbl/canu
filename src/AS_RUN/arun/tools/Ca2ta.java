// A quick hack of ca2ta.pl to output a .contig file given a .frg and .asm file

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.util.zip.*;

class StringCompression {
	public static final byte[] compress(String str) throws IOException {
		ByteArrayOutputStream out = new ByteArrayOutputStream();
		ZipOutputStream zout = new ZipOutputStream(out);
		zout.putNextEntry(new ZipEntry("0"));
		zout.write(str.getBytes());
		zout.closeEntry();
		byte[] compressed = out.toByteArray();
		zout.close();
		return compressed;
	}

	public static final String decompress(byte[] compressed) throws IOException {
		ByteArrayOutputStream out = new ByteArrayOutputStream();
		ByteArrayInputStream in = new ByteArrayInputStream(compressed);
		ZipInputStream zin = new ZipInputStream(in);
		zin.getNextEntry();
		byte[] buffer = new byte[1024];
		int offset = -1;
		while((offset = zin.read(buffer)) != -1) {
			out.write(buffer, 0, offset);
		}
		String decompressed = out.toString();
		out.close();
		zin.close();
		return decompressed;
	}

}

class Contig implements Serializable {
	String acc;
	int len;
	int npc;
	byte [] cns;
	ArrayList<MPS> mpsList;
	int [] ra_offsets;
	
	public String toString() {
	   	StringBuilder result = new StringBuilder();
		try {
			Comparator<MPS> comp = new Comparator<MPS>() { 
				public int compare(MPS o1, MPS o2) { 
					int retVal = o1.lpos - o2.lpos;
					if ( retVal == 0 ) {
						FRG frg1 = Ca2ta.frgHash.get(o1.mid);
						FRG frg2 = Ca2ta.frgHash.get(o2.mid);
						int o1_len = o1.dln+frg1.rclr-frg1.lclr+1;
						int o2_len = o2.dln+frg2.rclr-frg2.lclr+1;
						String o1_nm = frg1.nm;
						String o2_nm = frg1.nm;
						retVal = o1_len - o2_len;
					}
					return retVal;
				}};
			Collections.sort(mpsList,comp);
			String cnsDecompressed = StringCompression.decompress(cns);
			result.append("##" + acc + " " + npc + " " + cnsDecompressed.length() + " bases, 00000000 checksum." + "\n");
			result.append(Ca2ta.formatSeq(cnsDecompressed));
			buildOffsets();			
			for ( MPS mps : mpsList )
				result.append(printMPS(mps));
		} catch (java.io.IOException e) {
			e.printStackTrace();
		}
		return result.toString();
	}
	public void buildOffsets() throws IOException {
		String cnsStr = StringCompression.decompress(cns);
		ra_offsets = new int[cnsStr.length()];		
		int coord = 0;
		for (int i = 0; i < cnsStr.length(); i++){
			if ( cnsStr.charAt(i) != '-'){
				coord++;
			}
			ra_offsets[i] = coord;
		}
	}
	public String printMPS ( MPS mps ) throws IOException {
		StringBuilder result = new StringBuilder();

		FRG frg = Ca2ta.frgHash.get(mps.mid);
		if ( frg.seq == null ) { 
			return "";
		}
		String seqDecompressedFull = StringCompression.decompress(frg.seq);
		String seqDecompressed = seqDecompressedFull.substring(frg.lclr,frg.rclr);
		seqDecompressed = seqDecompressed.toUpperCase();
		frg.lclr++;
		int asml = mps.lpos; //asml
		int asmr = mps.rpos; //asmr
		int seqleft = frg.lclr; //seqleft
		int seqright = frg.rclr; //seqright

		if ( mps.rc ) {
			seqDecompressed = MPS.rc(seqDecompressed);
			seqleft = frg.rclr;
			seqright = frg.lclr;
		}

		int offset = asml;

		if ( mps.dln > 0 )
			seqDecompressed = mps.insertDeletes(seqDecompressed);

		if ( asmr - asml > seqDecompressed.length()) {
			asmr = asml + seqDecompressed.length() - 1;
		}
		if ( asmr <= 0 )
			asmr = ra_offsets[0];
		else
			asmr = ra_offsets[asmr - 1];
		if ( asml >= ra_offsets.length )
			asml = ra_offsets[ra_offsets.length-1];
		else
			asml = ra_offsets[asml];
		result.append("#" + frg.nm + "(" + offset + ")" 
			+ " [" + (mps.rc?"RC":"") + "] " 
			+ (mps.dln+frg.rclr-frg.lclr+1) 
			+ " bases, 00000000 checksum." 
			+ " {" + seqleft + " " + seqright + "}" 
			+ " <" + asml + " " + asmr +">" + "\n");
		result.append(Ca2ta.formatSeq(seqDecompressed));
		
		return result.toString();	
	}

}

class FRG implements Serializable {
	String mid;
	String nm;
	byte [] seq;
	int lclr;
	int rclr;
	
}

class MPS implements Serializable {
	boolean rc;
	char typ;
	String mid;
	int lpos;
	int rpos;
	int dln; //gapno
	int del[]; //gaps
	
	String insertDeletes(String seq) {
		StringBuffer retSB = new StringBuffer(seq);
		int localOffset = 0;
		for ( int gapPos : del ) {
			retSB.insert(gapPos+localOffset++,"-");
		}
		return retSB.toString();
	}
	
	static String rc(String seqStr) {
		StringBuffer retSB = new StringBuffer();

		for(int i=0 ; i<seqStr.length() ; i++) {
        		switch(seqStr.charAt(i)) {
        			case 'A': retSB.append('T');break;
                		case 'T': retSB.append('A');break;
                		case 'G': retSB.append('C');break;
                		case 'C': retSB.append('G');break;
                		default: retSB.append('N');break;
        		}
		}
		return retSB.reverse().toString();        
	}	
}

public class Ca2ta {
   static final int SEQ_OUTPUT_SIZE = 60;
   static Pattern acc = Pattern.compile("^acc:\\((\\d+)");
   static Pattern bracket = Pattern.compile("^\\{(\\w+)$");
   static Pattern field = Pattern.compile("^(\\w\\w\\w):(.*)$");

   static Map<String,FRG> frgHash;
   static Map<String,Contig> ccoHash;
   static String prefix;
   static String wrkDir;

   public static void main( String [] argv ) {
	if ( argv.length != 1 ) {
		System.err.println("Usage: Ca2ta <prefix>");
		System.err.println("       Requires a <prefix>.frg and <prefix>.asm files.");
		return;
	}
	File file = new File(argv[0]);
	prefix = file.getName();
	wrkDir = file.getParent();
	if ( wrkDir == null )
		wrkDir = ".";
	wrkDir = wrkDir.concat("/");
	System.out.println("prefix: " + prefix);
	System.out.println("wrkDir: " + wrkDir);

	if ( !frgHashExists() ) {
		//Step 0: Create clr file
		createCLR();

		//Step 1: Read prefix.clv file
		System.out.println("\nreadCLR");
		readClr();

		//Step 2: Read the FRG file
		System.out.println("\nreadFRG");
		readFrg();

		writeFrg();
	} else {
		readFrgHash();
	}
	
	if ( !ccoHashExists() ) {
		//Step 2: Read the noAFG file
		System.out.println("\nreadASM");
		readASM();	
		
		writeCCO();	
	} else {
		readCCOHash();
	}


	//Step 3: Print the contig file
	System.out.println("\nprintContig");
	printContig();

   }
   
   static void createCLR() {
   	//First check if it already exists
	String clrFile = wrkDir.concat(prefix).concat(".clr");
	
	if ( new File(clrFile).exists() )
		return;
   	
	try
        {   
		String cmd = "/usr/local/devel/ATG/moweis/CA_Latest/src/AS_RUN/arun/tools/asmToCLR.sh  " + wrkDir.concat(prefix).concat(".asm");
		System.out.println("Executing: " + cmd); //System.exit(1);
		Runtime rt = Runtime.getRuntime();
		Process proc = rt.exec(cmd);
		BufferedReader ls_in = new BufferedReader(
                        		  new InputStreamReader(proc.getInputStream()));
		FileWriter output = new FileWriter(clrFile);
		BufferedWriter bufWrite = new BufferedWriter(output);
		String ls_str;

		while ((ls_str = ls_in.readLine()) != null) {
		    output.write(ls_str+"\n");
		}
		ls_in.close();
		bufWrite.close();
        } catch (Throwable t) {
            t.printStackTrace();
        }     
	
   }
   
   static void readClr() {
	String line = null;    // String that holds current file line
	try {
		FileReader input = new FileReader(wrkDir.concat(prefix).concat(".clr"));

		BufferedReader bufRead = new BufferedReader(input);
		int count = 0;  // Line number of count 

		line = bufRead.readLine();

		frgHash = new HashMap<String,FRG>();

		while (line != null){
		  String[] values = line.split("\\s");
		  FRG frg = new FRG();
		  frg.mid = values[0];
		  frg.lclr = Integer.parseInt(values[1]);
		  frg.rclr = Integer.parseInt(values[2]);
		  
		  frgHash.put(frg.mid,frg);
		  line = bufRead.readLine();
		  count++;
		  if ( count % 100000 == 0 ) System.out.println("CLR: " + count);

		}
		System.out.println("CLR Final Read: " + count);
		System.out.println("CLR Hash size: " + frgHash.size());
		bufRead.close();              
	}catch (ArrayIndexOutOfBoundsException e){
		System.out.println("Error reading line: " + line + "	\n");
		e.printStackTrace();
	}catch (IOException e){
		e.printStackTrace();
	}
   }   

   static void readFrg() {
	try {
		FileReader input = new FileReader(wrkDir.concat(prefix).concat(".frg"));

		BufferedReader bufRead = new BufferedReader(input);
   		String line = bufRead.readLine();      

        	int count = 0;
		
		while ( line != null ) {
			Matcher m = bracket.matcher(line);
			
			if ( m.lookingAt() ) {
				String recName = m.group(1);
				if (recName.equals("FRG") ) { 
					FRG frg = processFRG(bufRead);
					if ( frg != null )  {
						FRG getFRG = frgHash.get(frg.mid);
						if ( getFRG != null ) {
							getFRG.nm = frg.nm;
							getFRG.seq = frg.seq;
							count++;
			   	        		if ( count % 100000 == 0 ) System.out.println("FRG: " + count);
						}
					} else {
						System.out.println("FRG IS NULL for: " + frg.mid);
					}
				} else {
					skipMsg(bufRead);
				}
			}
			line = bufRead.readLine();
		}

		System.out.println("FRG Final Read: " + count);
		System.out.println("FRG Hash size: " + frgHash.size());
		bufRead.close();              	
	}catch (ArrayIndexOutOfBoundsException e){
		System.out.println("Usage: java ReadFile filename\n");
		e.printStackTrace();
	}catch (IOException e){
		e.printStackTrace();
	}catch (Exception e){
		e.printStackTrace();
	}
   }

   static boolean frgHashExists() {
   	return new File(wrkDir.concat("frgHash.ser")).exists();
   }
   
   static boolean ccoHashExists() {
   	return new File(wrkDir.concat("ccoHash.ser")).exists();
   }

   @SuppressWarnings({"unchecked"})   
   static void readFrgHash () {
	//declared here only to ensure visibilty in finally clause
	ObjectInput input = null;
	try{
		//use buffering
		InputStream file = new FileInputStream( wrkDir.concat("frgHash.ser") );
		InputStream buffer = new BufferedInputStream( file );
		input = new ObjectInputStream ( buffer );
		//deserialize the List
		frgHash = (Map<String,FRG>)input.readObject();
	} catch(IOException e){
		e.printStackTrace();
	} catch (ClassNotFoundException e){
		e.printStackTrace();
	}
	finally{
		try {
			if ( input != null ) {
			  //close "input" and its underlying streams
			  input.close();
			}
		} catch (IOException e){
			e.printStackTrace();
		}
	}   
   }
   
   @SuppressWarnings({"unchecked"})   
   static void readCCOHash () {
	//declared here only to ensure visibilty in finally clause
	ObjectInput input = null;
	try{
		//use buffering
		InputStream file = new FileInputStream( wrkDir.concat("ccoHash.ser") );
		InputStream buffer = new BufferedInputStream( file );
		input = new ObjectInputStream ( buffer );
		//deserialize the List
		ccoHash = (Map<String,Contig>)input.readObject();
	} catch(IOException e){
		e.printStackTrace();
	} catch (ClassNotFoundException e){
		e.printStackTrace();
	}
	finally{
		try {
			if ( input != null ) {
			  //close "input" and its underlying streams
			  input.close();
			}
		} catch (IOException e){
			e.printStackTrace();
		}
	}   
   }
   
   static void writeFrg() {
	//declared here only to ensure visibilty in finally clause
	ObjectOutput output = null;
	try{
		//use buffering
		OutputStream file = new FileOutputStream( wrkDir.concat("frgHash.ser") );
		OutputStream buffer = new BufferedOutputStream( file );
		output = new ObjectOutputStream( buffer );
		output.writeObject(frgHash);
	} catch(IOException e){
		e.printStackTrace();
	} finally{
		try {
			if (output != null) {
			  //flush and close "output" and its underlying streams
			  output.close();
			}
		} catch (IOException e ){
			e.printStackTrace();
		}
	}   
   }

   static void writeCCO() {
	//declared here only to ensure visibilty in finally clause
	ObjectOutput output = null;
	try{
		//use buffering
		OutputStream file = new FileOutputStream( wrkDir.concat("ccoHash.ser") );
		OutputStream buffer = new BufferedOutputStream( file );
		output = new ObjectOutputStream( buffer );
		output.writeObject(ccoHash);
	} catch(IOException e){
		e.printStackTrace();
	} finally{
		try {
			if (output != null) {
			  //flush and close "output" and its underlying streams
			  output.close();
			}
		} catch (IOException e ){
			e.printStackTrace();
		}
	}   
   }
   
   
   static void readASM() {
	try {
		FileReader input = new FileReader(wrkDir.concat(prefix).concat(".asm"));
		BufferedReader bufRead = new BufferedReader(input);
		ccoHash = new HashMap<String,Contig>();
		getCARecord(bufRead);	
		bufRead.close();              
		System.out.println("Final CCO hash size: " + ccoHash.size());
	}catch (ArrayIndexOutOfBoundsException e){
		System.out.println("Usage: java ReadFile filename\n");
		e.printStackTrace();
	}catch (IOException e){
		e.printStackTrace();
	}catch (Exception e){
		e.printStackTrace();
	}   
   }
   
   static void printContig() {
		try {
			FileWriter output = new FileWriter(wrkDir.concat(prefix).concat(".contig"));
			BufferedWriter bufWrite = new BufferedWriter(output);
			writeContigs(bufWrite);	
			bufWrite.close();              
		}catch (ArrayIndexOutOfBoundsException e){
			System.out.println("Usage: java ReadFile filename\n");
			e.printStackTrace();
		}catch (IOException e){
			e.printStackTrace();
		}catch (Exception e){
			e.printStackTrace();
		}
   }   
   
   static void getCARecord ( BufferedReader bufRead ) throws Exception {
   	String line = bufRead.readLine();
	int count = 0;
	while ( line != null ) {
		Matcher m = bracket.matcher(line);
		String recName = "";
		if ( m.lookingAt() ) {
			recName = m.group(1);
		}

		if (recName.equals("CCO") ) { 
			processContig(bufRead);
			count++;
			if ( count % 10000 == 0 ) System.out.println("CCO: " + count);
		} else {
			skipMsg(bufRead);
		}
		line = bufRead.readLine();
	}        
   }
   
   static FRG processFRG ( BufferedReader bufRead ) throws Exception {
   
   	String line = bufRead.readLine();
	FRG frg = new FRG();
	boolean accSet = false;

	while ( line != null && !line.equals("}")) {
		Matcher fieldMatch = field.matcher(line);
		String fieldName = "";
		String fieldValue;
		try {
		if ( fieldMatch.lookingAt() )  {
			fieldName = fieldMatch.group(1);
			fieldValue = fieldMatch.group(2);			
			if ( !accSet && fieldName.equals("acc") ) {
				frg.mid = fieldValue;
				accSet = true;
			} else if ( fieldName.equals("src") ) {
				frg.nm = bufRead.readLine();
			} else if ( fieldName.equals("seq") ) {			
				frg.seq = StringCompression.compress(readSequence(bufRead));
				if ( frg.seq == null ) System.out.println("Null seq for frg: " + frg);
			}
		}
		} catch (NumberFormatException e)  {
			System.out.println("NumberFormatException at line: " + line);
			e.printStackTrace();
		}
		line = bufRead.readLine();
	}
	
	if ( frg.nm != null && frg.nm.equals(".") )
		frg.nm = frg.mid;
	return frg;   
   }

   static void processContig( BufferedReader bufRead ) throws Exception {
   	String line = bufRead.readLine();

	Matcher accMatch = acc.matcher(line);
	accMatch.lookingAt();
	Contig cco = new Contig();
	cco.acc = accMatch.group(1);	
	cco.mpsList = new ArrayList<MPS>();
	while ( line != null && !line.equals("}")) {
		Matcher bracketMatch = bracket.matcher(line);
		Matcher fieldMatch = field.matcher(line);
		String recName = "";
		String fieldName = "";
		String fieldValue;
		if ( bracketMatch.lookingAt() )
			recName = bracketMatch.group(1);
			if ( recName.equals("VAR") || recName.equals("UPS")) {
				skipMsg(bufRead);
			}
			else if ( recName.equals("MPS") ) {
				MPS returnedMPS = readMPS(bufRead);
				if ( returnedMPS != null ) {
					cco.mpsList.add(returnedMPS);
				}
			}
		else if ( fieldMatch.lookingAt() )  {
			fieldName = fieldMatch.group(1);
			fieldValue = fieldMatch.group(2);			
			if ( fieldName.equals("len") ) {
				cco.len = Integer.parseInt(fieldValue);
			} else if ( fieldName.equals("cns") ) {
				cco.cns = StringCompression.compress(readSequence(bufRead));
			} else if ( fieldName.equals("npc") ) {			
				cco.npc = Integer.parseInt(fieldValue);
			}
		} 
		line = bufRead.readLine();
	}
	ccoHash.put(cco.acc,cco);
   }
   
   static MPS readMPS ( BufferedReader bufRead) throws Exception {
   	String line = bufRead.readLine();
	
   	MPS mps = new MPS();
	while ( line != null && !line.equals("}")) {
		Matcher fieldMatch = field.matcher(line);
		String fieldName = "";
		String fieldValue;
		if ( fieldMatch.lookingAt() )  {
			fieldName = fieldMatch.group(1);
			fieldValue = fieldMatch.group(2);			
			if ( fieldName.equals("typ") ) {
				if ( !fieldValue.equals("R") ) {
					return null;
				}
				mps.typ = fieldValue.charAt(0);
			} else if ( fieldName.equals("mid") ) {
				mps.mid = fieldValue;				
			} else if ( fieldName.equals("pos") ) {			
				String [] posArray = fieldValue.split(",");
				mps.lpos = Integer.parseInt(posArray[0]);
				mps.rpos = Integer.parseInt(posArray[1]);
				if ( mps.lpos > mps.rpos ) {
					int temp = mps.rpos;
					mps.rpos = mps.lpos;
					mps.lpos = temp;
					mps.rc = true;	
				}
				
			} else if ( fieldName.equals("dln") ) {
				mps.dln = Integer.parseInt(fieldValue);
			} else if ( fieldName.equals("del") ) {
				line = bufRead.readLine();
				if ( line.equals("}") )  {
					break;
				}
				String delStr = new String();
				while ( line != null && !line.equals("}") ) {
					delStr = delStr.concat(line).concat(" ");
					line = bufRead.readLine();
				}
				String [] delStrArray = delStr.split(" ");
				mps.del = new int[delStrArray.length];
				for( int i= 0; i < delStrArray.length ; i++ )
					mps.del[i] = Integer.parseInt(delStrArray[i]);
				if ( line.equals("}") )
					break;
			}
		} 
		line = bufRead.readLine();
	}
	return mps;
   }
   
   static void skipMsg ( BufferedReader bufRead ) throws Exception {
   	String line = bufRead.readLine();
	int bracketCount = 1;
	while ( line != null) {
		if ( line.equals("}") ) {
			bracketCount--;
		}
		else {
			Matcher bracketMatch = bracket.matcher(line);
			if ( bracketMatch.lookingAt() )
				bracketCount++;		
		}
		
		if ( bracketCount == 0)
			return;
		line = bufRead.readLine();
	}   
   }
   
   
   static String readSequence ( BufferedReader bufRead ) throws Exception {
   	String line = bufRead.readLine();
	String sequence = "";
	
	while ( line != null && !line.equals(".")) {
		sequence = sequence.concat(line);
		line = bufRead.readLine();
	}	
	
   	return sequence;
   }  
   
   static void writeContigs ( BufferedWriter bufWrite ) throws Exception {
	for (String contigId :  new TreeSet<String>(ccoHash.keySet()))
		bufWrite.write(ccoHash.get(contigId).toString());
   }
   
   static String formatSeq (String seq) {
	StringBuilder result = new StringBuilder();
		
   	if ( seq != null )
   		for (int i = 0 ; i < seq.length() ; i+= SEQ_OUTPUT_SIZE ) {
			int diff = seq.length() - i;
			if ( diff <= 0 ) break;
			int endpoint = ( diff < SEQ_OUTPUT_SIZE ) ? diff : SEQ_OUTPUT_SIZE;
			result.append(seq.substring(i,i+endpoint) + "\n");
		}
	return result.toString();
   }
}
