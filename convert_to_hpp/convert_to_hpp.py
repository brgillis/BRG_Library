import sys
import os.path

def main(argv):
    # Read in command-line arguments
    if(len(argv)) <= 1:
        raise Exception("Name of data file must be passed at command-line.\n" + \
                    "eg. python convert_to_hpp.py my_table.dat")
        
    data_file_name = argv[1]
    name_root = os.path.splitext(data_file_name)[0]
    name = name_root[name_root.rfind('/')+1:-1]
    upper_name = name.upper()
    hpp_name = name_root + ".hpp"
    
    # Start writing
    with open(hpp_name,'w') as fo:
        # Write comment
        fo.write("/**\n * This file was automatically generated with convert_to_hpp.py, run on\n"
                 + " * " + data_file_name + "\n*/\n\n")
        
        # Write include guard (top)
        fo.write("#ifndef _BRG_" + upper_name + "_HPP_INCLUDED_\n" +
                 "#define _BRG_" + upper_name + "_HPP_INCLUDED_\n\n")
        
        # Write beginning of data storage
        fo.write("const char *"+name+"_data[] = {\n")
                                    
        # Write the data storage line by line
        with open(data_file_name,'r') as fi:
            for line in fi:
                fo.write("\""+line[0:-1] + "\",\n")

        # Write end of data storage
        fo.write("\"\\0\" };\n\n")
        
        # Write include guard (bottom)
        fo.write("#endif\n")
    


if __name__ == "__main__":
    main(sys.argv)