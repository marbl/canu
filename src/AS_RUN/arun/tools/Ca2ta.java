import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.lang.reflect.Field;

class ToString {
	/**
	* Intended only for debugging.
	*
	* <P>Here, a generic implementation uses reflection to print
	* names and values of all fields <em>declared in this class</em>. Note that
	* superclass fields are left out of this implementation.
	*
	* <p>The format of the presentation could be standardized by using
	* a MessageFormat object with a standard pattern.
	*/
	public String toString() {
	  StringBuilder result = new StringBuilder();
	  String newLine = System.getProperty("line.separator");

	  result.append( this.getClass().getName() );
	  result.append( "$ Object {" );
	  result.append(newLine);

	  //determine fields declared in this class only (no fields of superclass)
	  Field[] fields = this.getClass().getDeclaredFields();

	  //print field names paired with their values
	  for ( Field field : fields  ) {
	    result.append("  ");
	    try {
              result.append( field.getName() );
              result.append(": ");
              //requires access to private field:
              result.append( field.get(this) );
	    }
	    catch ( IllegalAccessException ex ) {
              System.out.println(ex);
	    }
	    result.append(newLine);
	  }
	  result.append("}");

	  return result.toString();
	}


}

class FRG extends ToString {
	String mid;
	String nm;
	String seq;
}

class Contig extends ToString {
	int len;
	int npc;
	String cns;
}

class MPS extends ToString {
	String typ;
	String mid;
	int lpos;
	int rpos;
	int dln; //gapno
	int del[]; //gaps
}
public class Ca2ta {
	static final int SEQ_OUTPUT_SIZE = 60;
	static Pattern acc = Pattern.compile("^acc:\\((\\d+)");
	static Pattern bracket = Pattern.compile("^\\{(\\w+)$");
	static Pattern field = Pattern.compile("^(\\w\\w\\w):(.*)$");
	
	static Map<String,String> clrHash;
	static Map<String,Contig> ccoHash;
	static Map<String,FRG> frgHash;
	static Map<String,ArrayList<MPS>> ccoMPSHash;
	static String prefix;

	public static void main( String [] argv ) {
		if ( argv.length != 1 ) return;
		else   prefix = argv[0];
		
		//Convert input .asm into .clv and .cco files		
		//Step 1: Read prefix.clv file
		System.out.println("readCLR\n");
		readClr();
			
		//Step 2: Read the FRG file
		System.out.println("readFRG\n");
		readFrg();
		
		
		//Step 2: Read the noAFG file
		System.out.println("readCCO\n");
		readCCO();

		//Step 3: Print the contig file
		System.out.println("printContig\n");
		printContig();

   }
   
   static void readFrg() {
	try {
		FileReader input = new FileReader(prefix.concat(".frg"));

		BufferedReader bufRead = new BufferedReader(input);
   		String line = bufRead.readLine();      

        	int count = 0;

		frgHash = new HashMap<String,FRG>();
		while ( line != null ) {
			Matcher m = bracket.matcher(line);
			String recName = "";
			if ( m.lookingAt() ) {
				recName = m.group(1);
			}
			int brackets = 0;
			if (recName.equals("FRG") ) { 
				FRG frg = processFRG(bufRead);
				if ( frg != null )  {
					frgHash.put(frg.mid,frg);
					count++;
			   	        if ( count % 100000 == 0 ) System.out.println("FRG: " + count);
				}
			} else {
				skipMsg(bufRead);
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
   
   static void readClr() {
	try {
		FileReader input = new FileReader(prefix.concat(".clr"));

		BufferedReader bufRead = new BufferedReader(input);
		String line;    // String that holds current file line
		int count = 0;  // Line number of count 

		line = bufRead.readLine();
		count++;

		clrHash = new HashMap<String,String>();

		while (line != null){
		  String[] values = line.split("\t");
		  clrHash.put(values[0],values[1].concat("\t").concat(values[2]));
		  line = bufRead.readLine();
		  count++;
		  if ( count % 100000 == 0 ) System.out.println("CLR: " + count);

		}
		System.out.println("CLR Final Read: " + count);
		System.out.println("CLR Hash size: " + clrHash.size());
		bufRead.close();              
	}catch (ArrayIndexOutOfBoundsException e){
		System.out.println("Usage: java ReadFile filename\n");
		e.printStackTrace();
	}catch (IOException e){
		e.printStackTrace();
	}
   }

   static void readCCO() {
		try {
			FileReader input = new FileReader(prefix.concat(".ccoSCF"));
			BufferedReader bufRead = new BufferedReader(input);
			ccoHash = new HashMap<String,Contig>();
			ccoMPSHash = new HashMap<String,ArrayList<MPS>>();
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
			FileWriter output = new FileWriter(prefix.concat(".contig.BETA"));
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

	while ( line != null ) {
		Matcher m = bracket.matcher(line);
		String recName = "";
		if ( m.lookingAt() ) {
			recName = m.group(1);
		}

		int brackets = 0;
		if (recName.equals("CCO") ) { 
			processContig(bufRead);
		} else {
			skipMsg(bufRead);
		}
		line = bufRead.readLine();
	}        
   }
   
   static FRG processFRG ( BufferedReader bufRead ) throws Exception {
   
   	String line = bufRead.readLine();
	FRG frg = new FRG();

	while ( line != null && !line.equals("}")) {
		Matcher fieldMatch = field.matcher(line);
		String fieldName = "";
		String fieldValue;
		if ( fieldMatch.lookingAt() )  {
			fieldName = fieldMatch.group(1);
			fieldValue = fieldMatch.group(2);			
			if ( fieldName.equals("acc") ) {
				frg.mid = fieldValue;
			} else if ( fieldName.equals("src") ) {
				frg.nm = bufRead.readLine();
			} else if ( fieldName.equals("seq") ) {			
				frg.seq = readSequence( bufRead );
			}
		} 
		line = bufRead.readLine();
	}
	return frg;   
   }
   
   static void processContig( BufferedReader bufRead ) throws Exception {
   	String line = bufRead.readLine();

	Matcher accMatch = acc.matcher(line);
	accMatch.lookingAt();
	String id = accMatch.group(1);
	Contig cco = new Contig();
	ArrayList<MPS> mpsList = new ArrayList<MPS>();
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
					mpsList.add(returnedMPS);
				}
			}
		else if ( fieldMatch.lookingAt() )  {
			fieldName = fieldMatch.group(1);
			fieldValue = fieldMatch.group(2);			
			if ( fieldName.equals("len") ) {
				cco.len = Integer.parseInt(fieldValue);
			} else if ( fieldName.equals("cns") ) {
				cco.cns = readSequence( bufRead );
			} else if ( fieldName.equals("npc") ) {			
				cco.npc = Integer.parseInt(fieldValue);
			}
		} 
		line = bufRead.readLine();
	}
	ccoHash.put(id,cco);   
	ccoMPSHash.put(id,mpsList);	
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
				mps.typ = fieldValue;
			} else if ( fieldName.equals("mid") ) {
				mps.mid = fieldValue;
			} else if ( fieldName.equals("pos") ) {			
				String [] posArray = fieldValue.split(",");
				mps.lpos = Integer.parseInt(posArray[0]);
				mps.rpos = Integer.parseInt(posArray[1]);
			} else if ( fieldName.equals("dln") ) {
				mps.dln = Integer.parseInt(fieldValue);
			} else if ( fieldName.equals("del") ) {
				line = bufRead.readLine();
				if ( line.equals("}") )  {
					return null;
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
   	Set<String> ccoSorted = new TreeSet<String>(ccoHash.keySet());
	for (String contigId :  ccoSorted) {
		Contig contig = ccoHash.get(contigId);
		bufWrite.write("##" + contigId + " " + contig.npc + " " + contig.cns.length() + " bases, 00000000 checksum.\n");
		printSeq(bufWrite,contig.cns);
		for ( MPS mps : ccoMPSHash.get(contigId) )
			bufWrite.write(mps.toString());
	}
   }
   
   static void printSeq (BufferedWriter bufWrite, String seq) throws Exception {
   	if ( seq == null ) return;
   	for (int i = 0 ; i < seq.length() ; i+= SEQ_OUTPUT_SIZE ) {
		int diff = seq.length() - i;
		if ( diff <= 0 ) return;
		int endpoint = ( diff < SEQ_OUTPUT_SIZE ) ? diff : SEQ_OUTPUT_SIZE;
		bufWrite.write(seq.substring(i,i+endpoint) + "\n");
	}   
   }
}
