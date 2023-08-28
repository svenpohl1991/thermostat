import pandas as pd
import numpy as np
import sys

def get_round_number(number,type):
    if type == 't':
        return round(float(number))
    elif type == 'p':
        if abs(float(number) - 0.101325) < 10E-12:
            return float(number)
        elif float(number) < 100.0:
            return round(float(number),1)
        elif float(number) >= 100.0:
            return round(float(number))

def rebuild_string_with_rounded_numbers(input_string,type):
    # Split the input string by the dash
    numbers = input_string.split('-')

    # Check if there are two numbers
    if len(numbers) == 2:
        try:
            # Convert the substrings to integers and round them
            num1 = get_round_number(numbers[0],type)
            num2 = get_round_number(numbers[1],type)

            # Rebuild the string with rounded numbers and return it
            result_string = f"{num1}-{num2}"
            return result_string
        except ValueError:
            return "Invalid input: Not numeric values"
    elif len(numbers) == 1:
        try:
            if type == 'p':
                return f"{get_round_number(numbers[0],type)}"
            else:
                return numbers[0]
        except ValueError:
            return "Invalid input: Not a numeric value"
    else:
        return "Invalid input: Not exactly one or two numbers"

class CSVFileParser:
    def __init__(self, file_path):
        self.file_path = file_path
        self.dataframes = []
        self.props = []
        self.mix_name = ''
        
    def parse_csv(self):
        current_content = []

        with open(self.file_path, 'r') as file:
            for line in file:
                if line.startswith('@'):
                    l = line.replace('#','')
                    l = l.replace('@','')
                    self.mix_name = l.replace('\n','').strip()
                
                if line.startswith('#'):
                    l = line.replace('#','')
                    l = l.replace('\n','').strip()
                    self.props.append(l)
                    
                    if current_content:
                        df = pd.DataFrame([line.strip().split(';') for line in current_content[1:]])
                        df.columns = current_content[0].strip().split(';')  # Use the first line as the header
                        self.dataframes.append(df)
                        current_content = []
                elif not '@' in line:
                    current_content.append(line)

        if current_content:
            df = pd.DataFrame([line.strip().split(',') for line in current_content[1:]])
            df.columns = current_content[0].strip().split(',')  # Use the first line as the header
            self.dataframes.append(df)

    def get_original_dataframe(self):
        return pd.read_csv(self.file_path)


class table_writer():
    def __init__(self):
        self.header_dict = { "TBOUNDS": "\makecell{$T_{\\text{min}} - T_{\\text{max}}$ \\\ K}",
                "AUTHOR": "\makecell{Author \\\ }",
                "PTSCALC" : "\makecell{$N$ \\\ }",
                "PBOUNDS": "\makecell{$p_{\\text{min}} - p_{\\text{max}}$ \\\  MPa}",
                "XBOUNDS": "\makecell{$x_{\\text{min}} - x_{\\text{max}}$ \\\  \%}",
                "YBOUNDS": "\makecell{$y_{\\text{min}} - y_{\\text{max}}$ \\\  \%}",
                "AARDY" : "\makecell{$\\text{AARD}_{\\text{y}}$ \\\ \%}",
                "AARDX" : "\makecell{$\\text{AARD}_{\\text{x}}$ \\\ \%}",
                "AARD" : "\makecell{$\\text{AARD}$ \\\ \%}",
                "PTSCALCFAILED" : "\makecell{$N_\\text{f}$ \\\ }"}

    def color_if_even(self,s):
        return ["background-color: red" if len(val) > 0 else '' for val in s]

    def get_caption(self,prop):
        if 'PVT' in prop:
            return r""

    def write_tables(self,file_path):
        parser = CSVFileParser(file_path)
        parser.parse_csv()
        self.dataframes = parser.dataframes
    
        # caption

    
        for i,data in enumerate(self.dataframes):
            data.columns = data.columns.str.replace(' ', '')
            
            # round numbers
            data_copy = []
            for j,t in enumerate(data["TBOUNDS"]):
                new = rebuild_string_with_rounded_numbers(t,'t')
                if not 'Invalid' in new:    
                        if not new in data_copy: 
                            data.loc[j,'TBOUNDS'] = new
                            data_copy.append(new)
                        else:
                            data.loc[j,'TBOUNDS'] = " "
                else:
                    data_copy = []

            for j,t in enumerate(data["PBOUNDS"]):
                new = rebuild_string_with_rounded_numbers(t,'p')
                if not 'Invalid' in new:
                    data.loc[j,'PBOUNDS'] = new


            data = data.rename(columns=self.header_dict)
            data = data.style.format(na_rep='')
            data = data.hide(axis="index")
            data.to_latex(parser.mix_name+'_'+parser.props[i]+'.tex',column_format='lccccccc',environment='longtable',hrules=True)

if __name__ == "__main__":
    t = table_writer()
    t.write_tables(sys.argv[1])




