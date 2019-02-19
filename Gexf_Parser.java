import java.io.File;
import java.util.ArrayList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.w3c.dom.Node;
import org.w3c.dom.Element;

public class Gexf_Parser {
	
	public Graph read_File(String fname){
		System.out.println("Reading " + fname + " file...");
		Graph g = new Graph();
		try{
			File inputFile = new File(fname);
			DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
			DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
			Document doc = dBuilder.parse(inputFile);
			doc.getDocumentElement().normalize();
			
			NodeList nList = doc.getElementsByTagName("node");
			NodeList eList = doc.getElementsByTagName("edge");
			
			g.al_of_nodes = new ArrayList<MyNode>();
			g.al_of_edges = new ArrayList<MyEdge>();
			
			
			for(int temp = 0; temp < nList.getLength(); temp++)
			{
				Node nNode = nList.item(temp);
				if(nNode.getNodeType() == Node.ELEMENT_NODE)
				{
					Element eElement = (Element) nNode;
					MyNode n = new MyNode();
					n.Id = eElement.getAttribute("id");
					n.label = eElement.getAttribute("label");
					g.al_of_nodes.add(n);					
				}				
				
			}
			for(int temp = 0; temp < eList.getLength(); temp++)
			{
				Node eNode = eList.item(temp);
				if(eNode.getNodeType() == Node.ELEMENT_NODE)
				{
					Element eElement = (Element) eNode;
					MyEdge e = new MyEdge();
					e.Id = eElement.getAttribute("id");
					e.source = eElement.getAttribute("source");
					e.target = eElement.getAttribute("target");
					e.weight = Double.parseDouble(eElement.getAttribute("weight"));
					g.al_of_edges.add(e);					
				}				
			}		
			
		}
		catch(Exception e){
			e.printStackTrace();
		}
		return g;
	}
	

}
