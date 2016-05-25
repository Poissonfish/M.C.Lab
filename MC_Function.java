import org.rosuda.JRI.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;
import net.miginfocom.swing.MigLayout;
import java.awt.ScrollPane;

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
	int loft = 2;
	int gs = 1;
	
	String Path ="";
	JButton start = new JButton("START");
	JTextField wd = new JTextField(20);

	JTextField folder = new JTextField(5);
	JTextField pf = new JTextField(20);
	JTextField pr = new JTextField(20);
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
	JLabel PF = new JLabel("Forward Sequence");
	JLabel PR = new JLabel("Reverse Sequence");
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
	
	Rengine r = new Rengine(new String[]{"--no-save"}, false, new console());
	
	public void submit(){
		String desdir = wd.getText();
		String fold = folder.getText();
		String primer_f = pf.getText();
		String primer_r = pr.getText();
		String name_primer_f = npf.getText();
		String name_primer_r = npr.getText();
		String local_path = lp.getText();
		String source = sc.getText();
		String uss = us.getText();
		String pdd = pw.getText();
		
		r.eval("folder='"+fold+"'");
		r.eval("primer_f='"+primer_f+"'");
		r.eval("primer_r='"+primer_r+"'");
		r.eval("name_primer_f='"+name_primer_f+"'");
		r.eval("name_primer_r='"+name_primer_r+"'");
		r.eval("local_path='"+local_path+"'");
		r.eval("source='"+source+"'");	
		r.eval("username='"+uss+"'");		
		r.eval("password='"+pdd+"'");
		r.eval("desdir='"+desdir+"'");

		if(loft==1){
			r.eval("local=FALSE");
		}else if (loft==2){
			r.eval("local=TRUE");
		}
		if(gs==-1){
			r.eval("nt_search=TRUE");
		}else if (gs==1){
			r.eval("nt_search=FALSE");
		}
		
		r.eval("setwd(desdir)");
		r.eval("source('./MCLab_Function.R')");
		r.eval("print(desdir)");
		r.eval("print(folder)");
		r.eval("time=proc.time()[3]");		
		r.eval("MCLab(primer_f, primer_r, name_primer_f, name_primer_r, source, username, password, desdir,folder,local_path, local, nt_search)");
		r.eval("proc.time()[3]-time");		
	}	

	public FrameLayout(){	
		wd.setText("/home/mclab/R/git/M.C.Lab");
		sc.setText("ftp://140.109.56.5/");
		us.setText("rm208");
		pw.setText("167cm");
	
		
		start.setPreferredSize(new Dimension(200,300));
		start.setFont(new Font("Arial", Font.BOLD, 40));
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
		
		JPanel prPanel= new JPanel(new MigLayout());
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
				"NCBI BLAST"));

		JPanel mainPanel= new JPanel(new MigLayout());
		mainPanel.setPreferredSize(new Dimension(590,480));
	
		mainPanel.add(srPanel,"cell 0 0 2 1, w 570::, dock north");		
		mainPanel.add(wdPanel,"cell 0 1, w 360::");		
		mainPanel.add(prPanel,"cell 0 2, w 360::");
		mainPanel.add(ruPanel,"cell 1 1 1 3, w 220::");			
		mainPanel.add(gePanel,"cell 0 3, h 70::,w 360::");
		
		
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
	    	gs=gs*(-1);  
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

class console implements RMainLoopCallbacks
{
    JFrame f;
	
    public JTextArea textarea = new JTextArea();	
    JScrollPane logScroll = new JScrollPane(textarea,ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS,ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
    
    public console() {   	
        f = new JFrame();
        f.setLocation(680,-100);
        f.getContentPane().add(new JScrollPane(textarea, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_NEVER));
        f.setSize(new Dimension(450,450));
       
        f.show();
    }
    
 	
    
    public void rWriteConsole(Rengine re, String text, int oType) {
        textarea.append(text);
    }
    
    public void rBusy(Rengine re, int which) {
        System.out.println("rBusy("+which+")");
    }
    
    public String rReadConsole(Rengine re, String prompt, int addToHistory) {
        System.out.print(prompt);
        try {
            BufferedReader br=new BufferedReader(new InputStreamReader(System.in));
            String s=br.readLine();
            return (s==null||s.length()==0)?s:s+"\n";
        } catch (Exception e) {
            System.out.println("jriReadConsole exception: "+e.getMessage());
        }
        return null;
    }
    
    public void rShowMessage(Rengine re, String message) {
        System.out.println("rShowMessage \""+message+"\"");
    }
    
    public String rChooseFile(Rengine re, int newFile) {
	FileDialog fd = new FileDialog(f, (newFile==0)?"Select a file":"Select a new file", (newFile==0)?FileDialog.LOAD:FileDialog.SAVE);
	fd.show();
	String res=null;
	if (fd.getDirectory()!=null) res=fd.getDirectory();
	if (fd.getFile()!=null) res=(res==null)?fd.getFile():(res+fd.getFile());
	return res;
    }
    
    public void   rFlushConsole (Rengine re) {
	}
    
    public void   rLoadHistory  (Rengine re, String filename) {
    }			
    
    public void   rSaveHistory  (Rengine re, String filename) {
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
  


