import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class Final_Score {
	public ArrayList<Score> read_final_scores(String fname) throws IOException{
		System.out.println("Reading " + fname + " file...");
		ArrayList<Score> al_of_score = new ArrayList<Score>();
		
		BufferedReader br = null;
		FileReader fr = null;
		
		try{
			fr = new FileReader(fname);
			br = new BufferedReader(fr);
			
			String Line;
			
			while((Line = br.readLine()) != null){
				Score s = new Score();
				s.n1 = new MyNode();
				s.n2 = new MyNode();
				
				String[] str_arr = Line.split(" ");
				s.n1.Id = str_arr[0];
				s.n1.label = str_arr[0];
				
				s.n2.Id = str_arr[1];
				s.n2.label = str_arr[1];
				
				s.score = Double.parseDouble(str_arr[2]);
				
				al_of_score.add(s);				
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
		return al_of_score;		
	}

}
