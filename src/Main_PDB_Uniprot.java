import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

public class Main_PDB_Uniprot {

	public Main_PDB_Uniprot() throws IOException {

		PrintStream csvWriter = null;
//		LinkedList<CsDatabaseEntry> allunidatanotcurated = new LinkedList<CsDatabaseEntry>();
		java.util.Calendar calendar = java.util.Calendar.getInstance();

		String version = "JUNE2012";
		int k = 0;
		int kFirst = 0;
		int kInt = 0;
	
		CsDatabaseEntry[] firstcapacityarray = null;
		CsDatabaseEntry[] intermediatecapacityarray = new CsDatabaseEntry[k];


		String UniprotIDListURL = "http://www.uniprot.org/uniprot/?query=(organism%3A10116+OR+organism%3A9606+OR+organism%3A10090)+reviewed%3Ayes%22+cleavage%3B%22&format=list";
		int lineCount = getLineCount(UniprotIDListURL);
		String UniprotURLdebut = "http://www.uniprot.org/uniprot/?query=(organism%3A10116+OR+organism%3A9606+OR+organism%3A10090)+reviewed%3Ayes%22+cleavage%3B%22&format=xml&limit=100&offset=";
//		String UniprotURL = "http://www.uniprot.org/uniprot/?query=(organism%3A10116+OR+organism%3A9606+OR+organism%3A10090)+reviewed%3Ayes%22+cleavage%3B%22&format=xml&limit=100";
//
		String UniprotURL = null;
		int offset = 0;
		while (offset < lineCount) {
			UniprotURL = UniprotURLdebut + offset;
			System.out.println(UniprotURL);
//			 String UniprotURL = "http://www.uniprot.org/uniprot/P01375.xml";

			NodeList entries = getEntries(
					"/uniprot/entry[./feature[@type='site'][contains(@description, 'Cleavage')]]",
					parseUniprot(UniprotURL));
			
			System.out.println(entries.getLength());
			for (int i = 0; i < entries.getLength(); i++) {

				SubstrateDatabaseEntry substratedatabase = new SubstrateDatabaseEntry();
				substratedatabase = getAccession(entries, i, substratedatabase);
				substratedatabase = getUniSubstratepproteinname(entries, i,
						substratedatabase);
				substratedatabase = getUniSubstrategenename(entries, i,
						substratedatabase);
				substratedatabase = getUniSubstratetaxonomy(entries, i,
						substratedatabase);
				substratedatabase = getUniSubstratesequence(entries, i,
						substratedatabase);

				String uniprotid = substratedatabase.getS_UniprotID();
				String genename = substratedatabase.getS_Symbol();
				String substrateTaxon = substratedatabase.getS_Taxon();
				String sequence = substratedatabase.getS_Sequence();
				System.out.println(sequence);
				String commentS = "-";

				NodeList csdbentries = getCsdbentries(
						"./feature[@type='site'][contains(@description, 'Cleavage')]",
						entries.item(i));
				

				if (csdbentries.getLength() > 0) {
					kFirst = csdbentries.getLength() + k ;
				} else if (csdbentries.getLength() == 0) {
					kFirst = k ;
				}
				firstcapacityarray = new CsDatabaseEntry[kFirst];
				System.arraycopy(intermediatecapacityarray, 0, firstcapacityarray, 0, k);
				System.out.println(csdbentries.getLength());
				int difference = 0;
				int l = 0;
				
				for (int j = 0; j < csdbentries.getLength(); j++) {
					
					firstcapacityarray[k+difference] = new CsDatabaseEntry();
					firstcapacityarray[k+difference].setSubstrate(substratedatabase);				
					
					Node n = csdbentries.item(j);
					
					String csdburl = "http://www.uniprot.org/uniprot/" + uniprotid;
					firstcapacityarray[k+difference].setExternal_Link(csdburl);
					SimpleDateFormat format = new SimpleDateFormat("dd-MM-yyyy");
	                Calendar originalDate = Calendar.getInstance();
	                String dateString = format.format(originalDate.getTime());
	                firstcapacityarray[k+difference].setCreation_Date(dateString);
	                firstcapacityarray[k+difference].setPMID("");
					
					
					LinkedList<String> p1intlist = getInformation("./location/begin/@position", n);
					
						
		                String p1int = null;
		                int intP1 = 0;
		                if (p1intlist.isEmpty()) {
		                    LinkedList<String> positionlist = getInformation("./location//position/@position", n);
		                    p1int = positionlist.getFirst();
		                    System.out.println(p1int);
		                    intP1 = Integer.parseInt(p1int);
		                    firstcapacityarray[k+difference].setP1_Position(intP1);
		                    char aaP1 = sequence.charAt(intP1 - 1);
		                    String saaP1 = Character.toString(aaP1);
		                    System.out.println(saaP1);
		                    firstcapacityarray[k+difference].setP1_Sequence(saaP1);
		                } else {
		                    p1int = p1intlist.getFirst();
		                    System.out.println(p1int);
		                    intP1 = Integer.parseInt(p1int);
		                    firstcapacityarray[k+difference].setP1_Position(intP1);
		                    char aaP1 = sequence.charAt(intP1 - 1);
		                    String saaP1 = Character.toString(aaP1);
		                    System.out.println(aaP1);
		                    firstcapacityarray[k+difference].setP1_Sequence(saaP1);
		                }

		                LinkedList<String> p1primeintlist = getInformation("./location/end/@position", n);
		                String p1primeint = null;
		                int intP1prime = 0;
		                if (p1primeintlist.isEmpty()) {
		                    intP1prime = intP1 + 1;
		                    System.out.println(intP1prime);
		                    firstcapacityarray[k+difference].setP1prime_Position(intP1prime);
		                    char aaP1prime = sequence.charAt(intP1prime - 1);
		                    String saaP1prime = Character.toString(aaP1prime);
		                    System.out.println(aaP1prime);
		                    firstcapacityarray[k+difference].setP1prime_Sequence(saaP1prime);

		                } else {
		                    p1primeint = p1primeintlist.getFirst();
		                    System.out.println(p1primeint);
		                    intP1prime = Integer.parseInt(p1primeint);
		                    firstcapacityarray[k+difference].setP1prime_Position(intP1prime);
		                    char aaP1prime = sequence.charAt(intP1prime - 1);
		                    String saaP1prime = Character.toString(aaP1prime);
		                    System.out.println(aaP1prime);
		                    firstcapacityarray[k+difference].setP1prime_Sequence(saaP1prime);
		                }

		                int sequencelength = sequence.length();
		                String cleavagesite = null;

		                if (intP1 - 3 > 0 || intP1 - 3 == 0) {
		                    if (intP1prime + 2 == sequencelength || intP1prime + 2 < sequencelength) {
		                        cleavagesite = sequence.substring(intP1 - 3, intP1prime + 2);
		                    } else {
		                    	if (intP1prime + 1 == sequencelength) {
		                        cleavagesite = sequence.substring(intP1 - 3, sequencelength) + "-";		                     
		                    	}
		                    	if (intP1prime == sequencelength) {
		                    	cleavagesite = sequence.substring(intP1 - 3, sequencelength) +"--";		                      
		                    	}
		                    }
		                } else {
		                    if (intP1prime + 2 == sequencelength || intP1prime + 2 < sequencelength) {
		                    	if (intP1 - 2 == 0) {
		                        cleavagesite = "-" + sequence.substring(0, intP1prime + 2);		                       
		                    	}
		                    	if (intP1 - 1 == 0) {
		                    	cleavagesite = "--" + sequence.substring(0, intP1prime + 2);			                
		                    	}
		                    } else {
		                        cleavagesite = sequence.substring(0, sequencelength);
		                        
		                    }
		                }
		                firstcapacityarray[k+difference].setCleavagesiteseaquence(cleavagesite);
		                System.out.println(cleavagesite);
					
					LinkedList<String> descriptionlist = getInformation(
							"./@description", n);
					String description = descriptionlist.getFirst();
					description = description.replaceAll("Cleavage; by ", "");
					description = description.replaceAll(
							"Cleavage, first; by ", "");
					description = description.replaceAll(
							"Cleavage, second; by ", "");
					description = description.replaceAll(",", ";");
					description = description.replaceAll("autolysis", genename);
					description = description.replaceAll("autocatalysis",
							genename);

					System.out.println(description);
					String proteaseTaxon = substrateTaxon;
//					proteasedatabase.setP_Taxon(proteaseTaxon);
					String curationUni = "-";
					
//					for (CsDatabaseEntry string : firstcapacityarray) {
//						System.out.println("MAIN " + string);
//					}
					
					intermediatecapacityarray = mapProteasetoLibrairy(commentS,
							proteaseTaxon, description, firstcapacityarray, k+difference, kFirst, curationUni);
					
					kInt = intermediatecapacityarray.length;
					l = kInt - firstcapacityarray.length ;		
					System.out.println("L " + l);
					if (l == 0) {
						k++;
					} else if ( l == 1) {
						difference++;
						difference++;
					} else if ( l == 2) {
						difference++;
						difference++;
						difference++;
					}
//					System.out.println("DIFFERENCE " + difference);
//					System.out.println("K " + k);
					firstcapacityarray = new CsDatabaseEntry[kInt];
					System.arraycopy(intermediatecapacityarray, 0, firstcapacityarray, 0, kInt);
					firstcapacityarray = intermediatecapacityarray;
//					for (CsDatabaseEntry string : firstcapacityarray) {
//						System.out.println("AFTER MAIN" + string);
//					}
//					firstcapacityarray = intermediatecapacityarray;
					System.out.println("NEW FIRST " + firstcapacityarray.length);
				}
				k = kInt;
			}
			
          offset = offset + 100;

			}

        try {
            System.out.println("-----------------");
            csvWriter = new PrintStream("UninotcuratedProteasixDB" + "_" + version + ".csv");
//            populateHeaders(csvWriter);
           for (CsDatabaseEntry csDatabaseEntry : intermediatecapacityarray) {
        	   System.out.println(csDatabaseEntry.getExternal_Link());
        	   System.out.println(csDatabaseEntry.getP1_Position());
        	   System.out.println(csDatabaseEntry.protease.getP_Symbol());
        	   System.out.println(csDatabaseEntry);
                populateData(csvWriter, csDatabaseEntry);
            }

        } catch (FileNotFoundException ex) {
            Logger.getLogger(Main_PDB_Uniprot.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            csvWriter.close();
        }
		
		}

	

	public static int getLineCount(String url) throws IOException {
		int linecounter = 0;
		BufferedReader urlReader = new BufferedReader(new InputStreamReader(
				new URL(url).openStream()));
		while (urlReader.readLine() != null) {
			linecounter++;
		}
		System.out.println(linecounter);
		return linecounter;
	}

	private NodeList getEntries(String query, Document xml) {
		XPathUniprot XPather = new XPathUniprot();
		NodeList entrylist = XPather.getNodeListByXPath(query, xml);
		return entrylist;
	}

	private Document parseUniprot(String url) {
		ParseUniprot parser = new ParseUniprot();
		Document xml = parser.getXML(url);
		xml.getXmlVersion();
		return xml;
	}

	private LinkedList<String> getInformation(String query, Node i) {
		XPathNodeUniprot XPathNoder = new XPathNodeUniprot();
		NodeList entrynodelist = XPathNoder.getNodeListByXPath(query, i);
		Loop l1 = new Loop();
		LinkedList<String> information = l1
				.getStringfromNodelist(entrynodelist);
		return information;
	}

	private SubstrateDatabaseEntry getAccession(NodeList entries, int i,
			SubstrateDatabaseEntry substratedatabase) {
		// GET URL AND GET ACCESSION using getInformation method
		LinkedList<String> idlist = getInformation("./accession[1]/text()",
				entries.item(i));
		String uniprotid = idlist.getFirst();
		System.out.println(uniprotid);
		substratedatabase.setS_UniprotID(uniprotid);
		return substratedatabase;
	}

	private SubstrateDatabaseEntry getUniSubstratepproteinname(
			NodeList entries, int i, SubstrateDatabaseEntry substratedatabase) {
		// GET SUBSTRATE PROTEIN NAME using getInformation method
		LinkedList<String> protnamelist = getInformation(
				"./protein/recommendedName/fullName/text()", entries.item(i));
		String protname = null;
		if (!protnamelist.isEmpty()) {
			protname = protnamelist.getFirst();
			protname = protname.replaceAll(",", "");
			System.out.println(protname);
			substratedatabase.setS_NL_Name(protname);
		}
		return substratedatabase;
	}

	private SubstrateDatabaseEntry getUniSubstrategenename(NodeList entries,
			int i, SubstrateDatabaseEntry substratedatabase) {
		// GET SUBSTRATE GENE NAME using getInformation method
		LinkedList<String> genenamelist = getInformation(
				"./gene/name[@type][1]/text()", entries.item(i));
		String genename = null;
		if (!genenamelist.isEmpty()) {
			genename = genenamelist.getFirst();
			System.out.println(genename);
			substratedatabase.setS_Symbol(genename);
		}
		return substratedatabase;
	}

	private SubstrateDatabaseEntry getUniSubstratetaxonomy(NodeList entries,
			int i, SubstrateDatabaseEntry substratedatabase) {
		// GET TAXONOMY using getInformation method
		String taxon = getInformation("./organism/dbReference/@id",
				entries.item(i)).getFirst();
		taxon.trim();
		if (taxon.equalsIgnoreCase("9606")) {
			taxon = "Homo Sapiens";
		} else if (taxon.equalsIgnoreCase("10116")) {
			taxon = "Rattus Norvegicus";
		}
		if (taxon.equalsIgnoreCase("10090")) {
			taxon = "Mus Musculus";
		}
		System.out.println(taxon);
		substratedatabase.setS_Taxon(taxon);
		return substratedatabase;
	}

	private SubstrateDatabaseEntry getUniSubstratesequence(NodeList entries,
			int i, SubstrateDatabaseEntry substratedatabase) {
		// GET PROTSEQUENCE using getInformation method
		LinkedList<String> sequencelist = getInformation("./sequence/text()",
				entries.item(i));
		if (!sequencelist.isEmpty()) {
		String sequence = sequencelist.getFirst();
		sequence = sequence.replaceAll("\n", "");
		substratedatabase.setS_Sequence(sequence);
		System.out.println(sequence);
		}
		return substratedatabase;
	}

	private NodeList getCsdbentries(String query, Node i) {
		XPathNodeUniprot XPathNoder = new XPathNodeUniprot();
		NodeList csdbentries = XPathNoder.getNodeListByXPath(query, i);
		return csdbentries;
	}

	private CsDatabaseEntry[] mapProteasetoLibrairy(String commentS, String proteaseTaxon,
			String proteaseDescription,
			CsDatabaseEntry[] firstcapacityarray, int j, int kFirst, String curationuni)
			throws IOException {
		
		CsDatabaseEntry intermediatecapacityarray[] = null;
		
		
		String commentP = commentS
				+ "; Check Protease Symbol and Accession; add to Substrate Librairy";
		firstcapacityarray[j].setComment(commentP);
		
		if (proteaseTaxon.equalsIgnoreCase("Homo Sapiens")) {
			BufferedReader bReader = createBufferedreader("/Users/julieklein/Dropbox/ProteasiX/LIBRAIRIES/ProteaseHSALibrairy.txt");
			intermediatecapacityarray = getProteaseInformation(bReader, proteaseDescription,
					firstcapacityarray, j, kFirst, proteaseTaxon, commentS, curationuni);
		} else {
			if (proteaseTaxon.equalsIgnoreCase("Rattus Norvegicus")) {
				BufferedReader bReader = createBufferedreader("/Users/julieklein/Dropbox/ProteasiX/LIBRAIRIES/ProteaseRNOLibrairy.txt");
				intermediatecapacityarray = getProteaseInformation(bReader, proteaseDescription,
						firstcapacityarray, j, kFirst, proteaseTaxon, commentS,
						curationuni);
			} else {
				if (proteaseTaxon.equalsIgnoreCase("Mus Musculus")) {
					BufferedReader bReader = createBufferedreader("/Users/julieklein/Dropbox/ProteasiX/LIBRAIRIES/ProteaseMMULibrairy.txt");
					intermediatecapacityarray = getProteaseInformation(bReader,
							proteaseDescription, firstcapacityarray, j, kFirst,
							proteaseTaxon, commentS, curationuni);
				}
			}
		}
		return intermediatecapacityarray;
	}

	private BufferedReader createBufferedreader(String datafilename)
			throws FileNotFoundException {
		BufferedReader bReader = new BufferedReader(
				new FileReader(datafilename));
		return bReader;

	}

	private CsDatabaseEntry[] getProteaseInformation(BufferedReader bReader,
			String proteaseDescription,
			CsDatabaseEntry[] firstcapacityarrray, int j, int kFirst,
			String proteaseTaxon, String commentS,
			String curationuni) throws IOException {
		String line;
		String commentP = null;
		int kInt = 0;
		String genename = null;
		String protname = null;
		
		Map<String, List<Set<String>>> hmap = new HashMap<String, List<Set<String>>>();
		CsDatabaseEntry intermediatecapacityarray[] = null;
		
		while ((line = bReader.readLine()) != null) {
			String splitarray[] = line.split("\t");
			String naturallanguage = splitarray[0];
			String uniId = splitarray[2];
			String ecnumber = splitarray[3];
			naturallanguage = naturallanguage.replaceAll("\"", "");
			String key = naturallanguage;
			if (!hmap.containsKey(key)) {
				List value = new ArrayList<Set<String>>();
				for (int j1 = 0; j1 < 2; j1++) {
					value.add(new HashSet<String>());
				}
				hmap.put(key, value);
			}
			hmap.get(key).get(0).add(naturallanguage);
			hmap.get(key).get(1).add(uniId + ";" + ecnumber);
		}
		
		ProteaseDatabaseEntry proteasedatabase = new ProteaseDatabaseEntry();
		proteasedatabase.setP_Taxon(proteaseTaxon);
		proteasedatabase.setP_NL_Name(proteaseDescription);
		proteasedatabase.setP_EC_Number("To check");
		proteasedatabase.setP_Name("To check");
		proteasedatabase.setP_Symbol("To check");
		proteasedatabase.setP_UniprotID("To check");
		firstcapacityarrray[j].setProtease(proteasedatabase);
		kInt = firstcapacityarrray.length;
		intermediatecapacityarray = new CsDatabaseEntry[kInt];
		System.arraycopy(firstcapacityarrray, 0, intermediatecapacityarray, 0, kFirst);
		intermediatecapacityarray = firstcapacityarrray;
		
		Iterator iterator = hmap.values().iterator();
		while (iterator.hasNext()) {
			String values = iterator.next().toString();
			String splitarray[] = values.split("\\], \\[");
			String naturallanguage = splitarray[0];
			naturallanguage = naturallanguage.replaceAll("\\[", "");
			if (naturallanguage.equals(proteaseDescription)) {
				String protease = splitarray[1].toString();
				String splitProtease[] = protease.split(", ");
				int splitPsize = splitProtease.length;
				if (splitPsize > 1) {
					kInt = firstcapacityarrray.length + splitPsize-1;
				} else {
					kInt = firstcapacityarrray.length;
				}
				System.out.println("j " + j);
				
				intermediatecapacityarray = new CsDatabaseEntry[kInt];
				System.arraycopy(firstcapacityarrray, 0, intermediatecapacityarray, 0, firstcapacityarrray.length);	
				ProteaseDatabaseEntry[] proteasedataarray = new ProteaseDatabaseEntry[splitPsize];
//				for (CsDatabaseEntry string : intermediatecapacityarray) {
//					System.out.println("SMALL " + string);
//				}
				for (int i = 0; i < splitPsize; i++) {		
					proteasedataarray[i] = new ProteaseDatabaseEntry();
					proteasedataarray[i].setP_Taxon(proteaseTaxon);
					proteasedataarray[i].setP_NL_Name(proteaseDescription);
					intermediatecapacityarray[j+i] = new CsDatabaseEntry();
					intermediatecapacityarray[j+i].setSubstrate(firstcapacityarrray[j].getSubstrate());
					intermediatecapacityarray[j+i].setExternal_Link(firstcapacityarrray[j].getExternal_Link());
					intermediatecapacityarray[j+i].setCreation_Date(firstcapacityarrray[j].getCreation_Date());
					intermediatecapacityarray[j+i].setPMID(firstcapacityarrray[j].getPMID());
					intermediatecapacityarray[j+i].setCleavagesiteseaquence(firstcapacityarrray[j].getCleavagesiteseaquence());
					intermediatecapacityarray[j+i].setP1_Position(firstcapacityarrray[j].getP1_Position());
					intermediatecapacityarray[j+i].setP1_Sequence(firstcapacityarrray[j].getP1_Sequence());
					intermediatecapacityarray[j+i].setP1prime_Position(firstcapacityarrray[j].getP1prime_Position());
					intermediatecapacityarray[j+i].setP1prime_Sequence(firstcapacityarrray[j].getP1prime_Sequence());
					intermediatecapacityarray[j+i].setCuration_Status(curationuni);
					String splitsplit[] = splitProtease[i].toString()
							.split(";");
					String uni = splitsplit[0].toString();
					System.out.println(uni + "uni");
					proteasedataarray[i].setP_UniprotID(uni);
					String EC = splitsplit[1].toString();
					EC = EC.replaceAll("\\]", "");
					proteasedataarray[i].setP_EC_Number(EC);
					
					if(!uni.equals("n.d.")) {
					String UniprotURL = "http://www.uniprot.org/uniprot/" + uni
							+ ".xml";
					NodeList entries = getEntries("/uniprot/entry",
							parseUniprot(UniprotURL));
					for (int j1 = 0; j1 < entries.getLength(); j1++) {
						 protname = getUniProteasepproteinname(entries, j1, proteasedataarray, i);
						 genename = getUniProteasegenename(entries, j1,
								 proteasedataarray, i);
						
					}
					proteasedataarray[i].setP_Name(protname);
					proteasedataarray[i].setP_Symbol(genename);
					intermediatecapacityarray[j+i].setProtease(proteasedataarray[i]);
					} else {
						proteasedataarray[i].setP_Name("n.d.");
						proteasedataarray[i].setP_Symbol("n.d.");
						intermediatecapacityarray[j+i].setProtease(proteasedataarray[i]);
					}
					commentP = commentS + ";-";
					intermediatecapacityarray[j+i].setComment(commentP);
				}
			}	
		}
//		for (CsDatabaseEntry csDatabaseEntry : intermediatecapacityarray) {
//			System.out.println("AFTER SMALL " + csDatabaseEntry);
//		}
		return intermediatecapacityarray;
	}

	public static void main(String[] args) throws MalformedURLException,
			IOException {
		// TODO code application logic here
		Main_PDB_Uniprot main = new Main_PDB_Uniprot();
	}

	private String getUniProteasepproteinname(NodeList entries, int i,
			ProteaseDatabaseEntry[] proteasedatabase, int j) {
		// GET SUBSTRATE PROTEIN NAME using getInformation method
		LinkedList<String> protnamelist = getInformation(
				"./protein/recommendedName/fullName/text()", entries.item(i));
		String protname = null;
		if (!protnamelist.isEmpty()) {
			protname = protnamelist.getFirst();
			protname = protname.replaceAll(",", "");
			System.out.println(protname);
//			proteasedatabase.setP_Name(protname);
		}
		return protname;
	}

	private String getUniProteasegenename(NodeList entries, int i,
			ProteaseDatabaseEntry[] proteasedatabase, int j) {
		// GET SUBSTRATE GENE NAME using getInformation method
		LinkedList<String> genenamelist = getInformation(
				"./gene/name[@type][1]/text()", entries.item(i));
		String genename = null;
		if (!genenamelist.isEmpty()) {
			genename = genenamelist.getFirst();
			System.out.println(genename);
//			proteasedatabase.setP_Symbol(genename);
		}
		return genename;
	}
	  private void populateData(PrintStream csvWriter, CsDatabaseEntry csdatabase) {
	        //System.out.println(cleavageSiteDBEntry);
	        
	        csvWriter.print(csdatabase.protease.getP_NL_Name());
	        csvWriter.print(",");
	        csvWriter.print(csdatabase.protease.getP_Symbol());
	        csvWriter.print(",");
	        csvWriter.print(csdatabase.protease.getP_UniprotID());
	        csvWriter.print(",");
	        csvWriter.print(csdatabase.protease.getP_EC_Number());
	        csvWriter.print(",");
	        csvWriter.print(csdatabase.protease.getP_Taxon());
	        csvWriter.print(",");
	        csvWriter.print(csdatabase.substrate.getS_NL_Name());
	        csvWriter.print(",");
	        csvWriter.print(csdatabase.substrate.getS_Symbol());
	        csvWriter.print(",");
	        csvWriter.print(csdatabase.substrate.getS_UniprotID());
//	                csvWriter.print(",");
//	                csvWriter.print(csdatabase.substrate.getSubstratesequence());
	        csvWriter.print(",");
	        csvWriter.print(csdatabase.substrate.getS_Taxon());
	        csvWriter.print(",");
	        csvWriter.print(csdatabase.getCleavagesiteseaquence());
	        csvWriter.print(",");
	        csvWriter.print(csdatabase.getP1_Position());
	        csvWriter.print(",");
	        csvWriter.print(csdatabase.getP1prime_Position());
	        csvWriter.print(",");
	        csvWriter.print(csdatabase.getP1_Sequence());
	        csvWriter.print(",");
	        csvWriter.print(csdatabase.getP1prime_Sequence());
	        csvWriter.print(",");
	        csvWriter.print(csdatabase.getExternal_Link());
	        csvWriter.print(",");
	        csvWriter.print(csdatabase.getPMID());
	        csvWriter.print(",");
	        csvWriter.print(csdatabase.getComment());
	        csvWriter.print(",");
	        csvWriter.print(csdatabase.getCuration_Status());
	        csvWriter.print(",");
	        csvWriter.print(csdatabase.getCreation_Date());
	        csvWriter.print("\n");
	    }

}