1. Install Python

Go to python.org

Download and install the newest version of Python 3

Important: During installation, check the box that says:

Add Python to PATH

Then click Install Now

2. Download the project from GitHub

Go to the GitHub page for the project

Click the green Code button

Click Download ZIP

3. Extract the ZIP file

Go to your Downloads folder

Right click the downloaded ZIP file

Click Extract All

After it finishes extracting, open the extracted project folder

You should see files like:

main.py
tree_builder.py
requirements.txt
birds.aln-clustalw

4. Open Command Prompt inside the project folder

Make sure you are inside the extracted project folder in File Explorer

At the top of File Explorer, click the address bar. This is the bar that shows the folder path, such as:

Downloads > Phylogenetic-Tree-Building-Analysis-main

Click inside that bar, delete the text and type:

cmd

A black Command Prompt window should open already inside the project folder.

To confirm you are in the right place, the command line should end with something like:

Phylogenetic-Tree-Building-Analysis-main>

5. Install the required packages

In Command Prompt, type:

python -m pip install -r requirements.txt

Wait for the installation to finish.

If the program later says a package is missing, run this command too:

python -m pip install biopython matplotlib pandas numpy seaborn cogent3 piqtree dendropy

6. Run the program

In the same Command Prompt window, type:

python main.py


7. Use these example inputs

When the program asks for the sequence alignment file, type:

birds.aln-clustalw

When it asks for the number of taxa or sequences, type:

4

When it asks for the sequence length, type:

282

When it asks whether you prefer speed or accuracy, type:

speed

The program should recommend UPGMA.

8. Build and visualize the tree

When it asks:

Would you like to build the recommended tree now?

Type:

yes

When it asks:

Would you like to visualize the tree?

Type:

yes

A window should open showing the generated phylogenetic tree.

9. How to run the program again later

Open the extracted project folder.

Click the address bar at the top of File Explorer.

Type:

cmd

Press enter

Then type:

python main.py