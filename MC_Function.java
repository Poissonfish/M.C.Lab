import java.util.Enumeration;
import org.rosuda.JRI.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.JFileChooser;
import javax.swing.border.*;
import net.miginfocom.swing.MigLayout;


public class MC_Function {
	public static void main(String[] args){        
		FrameLayout main = new FrameLayout();
		main.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		main.setTitle("M.C.Lab");		
		main.pack();
		main.show();
	}
}


class FrameLayout extends JFrame implements ActionListener{
	int loft = 0;
	int gs = 0;
	
	String Path ="";
	JButton start = new JButton("START");
	JTextField wd = new JTextField(20);
	JTextField folder = new JTextField(5);
	JTextField pf = new JTextField(10);
	JTextField pr = new JTextField(10);
	JTextField npf = new JTextField(3);
	JTextField npr = new JTextField(3);
	JTextField lp = new JTextField(20);
	JTextField us = new JTextField(10);
	JTextField pw = new JTextField(10);
	JTextField sc = new JTextField(20);
	JRadioButton buttononline = new JRadioButton("Online(FTP)");
	JRadioButton buttonlocal = new JRadioButton("Local");
	JCheckBox genesearch = new JCheckBox("BLAST sequences");
	JLabel WD = new JLabel("Working Directory");
	JLabel FD = new JLabel("Folder Name");
	JLabel PF = new JLabel("Forward Primer (Sequence)");
	JLabel PR = new JLabel("Reverse Primer (Sequence)");
	JLabel NPF = new JLabel("Forward Suffix");
	JLabel NPR = new JLabel("Reverse suffix");
	JLabel FP = new JLabel("Local File Path");
	JLabel US = new JLabel("Username");
	JLabel PW = new JLabel("Password");
	JLabel FT = new JLabel("Address URL");
	
	JButton browse = new JButton("Browse");
	JFileChooser choose = new JFileChooser();
	JButton browse2 = new JButton("Browse");
	JFileChooser choose2 = new JFileChooser();
	
	int value;
	int value2;
	
	Rengine r = new Rengine(new String[]{"--no-save"}, false, new TextConsole());
	
	public void submit(){
		// Path = central.getText();
		 r.eval("print('hello')");
		 System.out.println(loft);
		 System.out.println(gs);
	}
	

	public FrameLayout(){	
		buttonlocal.setSelected(true);
		us.setEnabled(false);
		pw.setEnabled(false);
		sc.setEnabled(false);
		ButtonGroup ftplocal = new ButtonGroup();
		ftplocal.add(buttononline);
		ftplocal.add(buttonlocal);
					
		start.addActionListener(this);
		browse.addActionListener(this);
		browse2.addActionListener(this);
		buttononline.addActionListener(this);
		buttonlocal.addActionListener(this);
		genesearch.addActionListener(this);		
		
		JPanel ruPanel = new JPanel();
		ruPanel.add(start);
		ruPanel.setBorder(new TitledBorder(new EtchedBorder(),
				"Run"));
		
		JPanel wdPanel = new JPanel(new MigLayout());
		wdPanel.add(WD,"cell 0 0");
		wdPanel.add(wd,"cell 0 1");
		wdPanel.add(browse,"cell 1 1");
		wdPanel.add(FD,"cell 0 2");
		wdPanel.add(folder,"cell 0 3");
		wdPanel.setBorder(new TitledBorder(new EtchedBorder(),
				"Destination"));	
		
		JPanel srPanel= new JPanel(new MigLayout());
		srPanel.add(buttononline, "cell 0 0");
		srPanel.add(FT, "cell 1 0");
		srPanel.add(sc, "cell 1 1");
		srPanel.add(US, "cell 2 0");
		srPanel.add(us, "cell 3 0");
		
		srPanel.add(PW, "cell 2 1");
		srPanel.add(pw, "cell 3 1");
		
		srPanel.add(buttonlocal, "cell 0 2");
		srPanel.add(FP, "cell 1 2");
		srPanel.add(lp, "cell 1 4");
		srPanel.add(browse2, "cell 2 4");
		srPanel.setBorder(new TitledBorder(new EtchedBorder(),
				"Data Source"));
		
		JPanel prPanel= new JPanel(new MigLayout("","[] 20 []",""));
		prPanel.add(PF, "cell 0 0");
		prPanel.add(pf, "cell 0 1");
		prPanel.add(NPF, "cell 1 0");
		prPanel.add(npf, "cell 1 1");
		prPanel.add(PR, "cell 0 2");
		prPanel.add(pr, "cell 0 3");
		prPanel.add(NPR, "cell 1 2");
		prPanel.add(npr, "cell 1 3");
		prPanel.setBorder(new TitledBorder(new EtchedBorder(),
				"Primers Information"));

		JPanel gePanel= new JPanel(new MigLayout());
		gePanel.add(genesearch);
		gePanel.setBorder(new TitledBorder(new EtchedBorder(),
				"BLAST"));

		JPanel mainPanel= new JPanel(new MigLayout());
		mainPanel.setPreferredSize(new Dimension(1200,600));
	
		mainPanel.add(wdPanel,"cell 0 0");
		mainPanel.add(ruPanel,"cell 1 0");
		mainPanel.add(srPanel,"cell 0 1");
		mainPanel.add(prPanel,"cell 0 2");
		mainPanel.add(gePanel,"cell 0 3");
	
		
		this.setContentPane(mainPanel);		
	}
	
	public void actionPerformed(ActionEvent ae){
	      Object source = ae.getSource();	
	      if (source == start){
	        submit();
	      }else if (source == buttononline){
	    	loft=1;  
	    	lp.setEnabled(false);
	    	us.setEnabled(true);
			pw.setEnabled(true);
			sc.setEnabled(true);
	      }else if (source == buttonlocal){
	        loft=2;
	        lp.setEnabled(true);
	    	us.setEnabled(false);
			pw.setEnabled(false);
			sc.setEnabled(false);
	      }else if (source == genesearch){
	    	gs=1;  
	      }else if (source == browse){
	    	choose.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
		    value= choose.showOpenDialog(null);
		    if (value == JFileChooser.APPROVE_OPTION){
		        File selectedfile = choose.getSelectedFile();
		    	wd.setText(selectedfile.getAbsolutePath());
		   	};
	      }else if (source == browse2){
		   	choose2.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
			value2= choose2.showOpenDialog(null);
			if (value2 == JFileChooser.APPROVE_OPTION){
				File selectedfile = choose2.getSelectedFile();
			    lp.setText(selectedfile.getAbsolutePath());
			};
		 }
	 }
}
	
/*
JFileChooser fileChooser = new JFileChooser();
//fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
int returnValue = fileChooser.showOpenDialog(null);
if (returnValue == JFileChooser.APPROVE_OPTION) {
  File selectedFile = fileChooser.getSelectedFile();
  System.out.println(selectedFile.getName());
  //File selectedFile = fileChooser.getCurrentDirectory();
  System.out.println(selectedFile.getAbsolutePath());
*/
  


