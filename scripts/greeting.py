import time


class Greeting:
    def __init__(self):
        # ASCII benzene molecule and MIHE information
        self.benzene_greet = '''
    **************************************************************************************
       /‾\   IIME: In silico Interface Mutagenesis Experiment tool
      | O |  Developed by M. Jukič, S. Kralj
       \_/   UM FKKT- Laboratory for Physical Chemistry and Chemical Thermodynamics
             UP FAMNIT- Jukič Laboratory for Cheminformatics
    **************************************************************************************   
    '''

    def start(self):
        # Print the greeting
        print(self.benzene_greet)
        time.sleep(3)

