import java.awt.BorderLayout;
import java.awt.Container;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.border.EtchedBorder;
import javax.swing.border.TitledBorder;

public class BroadLayout {

  public static void main(String[] args) {
    JFrame aWindow = new JFrame("This is a Box Layout");
    aWindow.setBounds(200, 200, 200, 200);
    aWindow.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

    Box left = Box.createVerticalBox();
    left.add(Box.createVerticalStrut(30));
    ButtonGroup radioGroup = new ButtonGroup();
    JRadioButton rbutton;
    radioGroup.add(rbutton = new JRadioButton("Red"));
    left.add(rbutton);
    left.add(Box.createVerticalStrut(30));
    radioGroup.add(rbutton = new JRadioButton("Green"));
    left.add(rbutton);
    left.add(Box.createVerticalStrut(30));
    radioGroup.add(rbutton = new JRadioButton("Blue"));
    left.add(rbutton);
    left.add(Box.createVerticalStrut(30));
    radioGroup.add(rbutton = new JRadioButton("Yellow"));
    left.add(rbutton);

    left.add(Box.createGlue());

    JPanel leftPanel = new JPanel(new BorderLayout());
    leftPanel.setBorder(new TitledBorder(new EtchedBorder(), "Line Color"));
    leftPanel.add(left, BorderLayout.CENTER);

    Box right = Box.createVerticalBox();
    right.add(Box.createVerticalStrut(30));
    right.add(new JCheckBox("Dashed"));
    right.add(Box.createVerticalStrut(30));
    right.add(new JCheckBox("Thick"));
    right.add(Box.createVerticalStrut(30));
    right.add(new JCheckBox("Rounded"));

    right.add(Box.createGlue());

    JPanel rightPanel = new JPanel(new BorderLayout());
    rightPanel
        .setBorder(new TitledBorder(new EtchedBorder(), "Line Properties"));
    rightPanel.add(right, BorderLayout.CENTER);

    Box top = Box.createHorizontalBox();
    top.add(leftPanel);
    top.add(Box.createHorizontalStrut(5));
    top.add(rightPanel);

    Container content = aWindow.getContentPane();
    content.setLayout(new BorderLayout());
    content.add(top, BorderLayout.CENTER);

    BoxLayout box = new BoxLayout(content, BoxLayout.Y_AXIS);

    content.setLayout(box);
    content.add(top);

    aWindow.setVisible(true);
  }
}