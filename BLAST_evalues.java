import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class BLAST_evalues {
	
	public ArrayList<Hit> read_Evalues(String fname) throws IOException{
		System.out.println("Reading " + fname + " file...");
		ArrayList<Hit> al_of_hit = new ArrayList<Hit>();
		
		BufferedReader br = null;
		FileReader fr = null;
		
		try{
			fr = new FileReader(fname);
			br = new BufferedReader(fr);
			
			String Line;
			
			while((Line = br.readLine()) != null){
				Hit h = new Hit();
				String[] str_arr = Line.split("\t");
				h.Node1 = str_arr[0];
				h.Node2 = str_arr[1];
				h.evalue = Double.parseDouble(str_arr[2]);
				
				al_of_hit.add(h);				
			}
			
		}
		catch(Exception e){
			e.printStackTrace();
		}
		finally{
			if(br != null)
				br.close();
			if(fr != null)
				fr.close();
		}		
		return al_of_hit;		
	}
}
class Hit{
	String Node1;
	String Node2;
	double evalue;
}
