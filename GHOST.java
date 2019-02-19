import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

public class GHOST {

	public static void main(String[] args) throws IOException {		
		
		Gexf_Parser gfp = new Gexf_Parser();

		Graph g1 = gfp.read_File("scere05.gexf");
		Graph g2 = gfp.read_File("scerehc.gexf");
		
		//Graph g1 = gfp.read_File("athal.gexf");
		//Graph g2 = gfp.read_File("dmel.gexf");
		
		//Graph g1 = gfp.read_File("cjejuni.gexf");
		//Graph g2 = gfp.read_File("ecoli.gexf");
		
		//Parameters
		int k = 2;
		double alpha = 1.0;
		String eval_file_name = "CJejuni_vs_EColi.evalues";
		double default_eval = 0.5;
		
		char controller = 'E';
		
		if(controller == 'S')
		{
			SeedPhase s1 = new SeedPhase();
			ArrayList<Score> al_of_s = s1.get_dist_among_all_pairs(g1, g2, k, alpha, eval_file_name, default_eval);
			Utility u = new Utility();
			u.write_score_file(al_of_s);
		}
		else if(controller == 'E')
		{
			Final_Score fs = new Final_Score();
			ArrayList<Score> al_of_s = fs.read_final_scores("Final_Scores.txt");
			Collections.sort(al_of_s, new Comparator<Score>() {
				@Override
		        public int compare(Score s1, Score s2)
		        {
		            return Double.compare(s1.score, s2.score);
		        }
			});
			
			Seed s = new Seed();
			s.score = 100000;
			for(Score score_elt : al_of_s)
			{
				if(score_elt.score < s.score)
				{
					s.n1 = score_elt.n1;
					s.n2 = score_elt.n2;
					s.score = score_elt.score;
				}
			}
			
			ExtendPhase e = new ExtendPhase();
			
			ArrayList<Alignment> al_of_a = e.extend_alignment(s, g1, g2, al_of_s, k, 0);

			Utility u = new Utility();
			u.write_alignment_file(al_of_a);		
			System.out.println("Total aligned node pairs: " + al_of_a.size());
		}
	}
}
