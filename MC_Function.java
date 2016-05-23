import java.util.Enumeration;
import org.rosuda.JRI.*;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

public class MC_Function extends JFrame{

	String Path ="";
	JTextField wd = new JTextField(10);
	JTextField folder = new JTextField(10);
	JTextField pf = new JTextField(10);
	JTextField pr = new JTextField(10);
	JTextField npf = new JTextField(3);
	JTextField npr = new JTextField(3);
	JTextField lp = new JTextField(10);
	JTextField us = new JTextField(10);
	JTextField pw = new JTextField(10);
	JTextField sc = new JTextField(10);
	
	JTextField central = new JTextField(10);

	Rengine r = new Rengine(new String[]{"--no-save"}, false, new TextConsole());
	
	public void submit(){
		 Path = central.getText();
		 r.eval(Path);
	}
	
	public MC_Function(){
		super("iPat");
		setSize(500,300);	
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		Container content = getContentPane();
		
		SpringLayout layout = new SpringLayout();
		content.setLayout(layout);	
		
		Component left = new JLabel("Specify your working directroy");	
		JButton right = new JButton("Run");
						
		content.add(right);
		content.add(central);
		content.add(left);
		right.addActionListener( (e)-> {
            submit();
        });			
		
		layout.putConstraint(SpringLayout.NORTH, left, 50, SpringLayout.NORTH, content);
		layout.putConstraint(SpringLayout.NORTH, central, 50, SpringLayout.NORTH, content);
		layout.putConstraint(SpringLayout.NORTH, right, 50, SpringLayout.NORTH, content);
		layout.putConstraint(SpringLayout.WEST, left, 10, SpringLayout.WEST, content);
		layout.putConstraint(SpringLayout.WEST, central, 10, SpringLayout.EAST, left);
		layout.putConstraint(SpringLayout.WEST, right, 10, SpringLayout.EAST, central);		
	}
		
	public static void main(String[] args){        
		//iPatwindow.setVisible(true); 
	    new MC_Function().setVisible(true);
	}
}	


/*
desdir="/home/mclab/R/git/M.C.Lab" 
setwd(desdir)
source('./MCLab_Function.R') 

time=proc.time()[3]
#######################
folder='ch11-11'
#11-11
primer_f="GGTGCTCAAGGCGGGACATTCGTT"
primer_r="TCGGGATTGGCCACAGCGTTGAC"
name_primer_f='T7'
name_primer_r='SP6'
#######################
local=TRUE
local_path='/home/mclab/R/git/M.C.Lab/seq/ch11-11/raw'
#######################
source="ftp://140.109.56.5/104DATA/0519/"
username="rm208"
password="167cm"
#######################
nt_search=TRUE
*/
