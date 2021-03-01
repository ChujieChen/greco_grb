class PATH:
    
    def __init__(self):
        try:
            import os
            import sys
            username = os.environ.get('USER')
        except:
            print("Cannot get USER right now.\n")
            raise

        icdata_dir = ""
        data_dir = ""
        ana_dir = ""

        if username == "cchen641":
            print("\n#### Working at GT ####\n")
            icdata_dir = "/storage/home/hhive1/cchen641/data/icecube/data/analyses/online_GRECO/version-002-p04"
            data_dir = "/storage/home/hhive1/cchen641/data/icecube/data/greco_grb/data"
            ana_dir = data_dir + "/csky_output"

        elif username == "cjchen":
            print("\n#### Running at IceCube ####\n")
            icdata_dir = "/data/user/cjchen/csky/analyses/online_GRECO/version-002-p04"
            data_dir = "/data/user/cjchen/2020-2021/Spring2021/greco_grb/data"
            ana_dir = data_dir + "/csky_output" 
        else: 
            print("#### USER not recognized ####\n")
        
        self.USER = username
        self.ICDATA_DIR = icdata_dir
        self.DATA_DIR = data_dir
        self.ANA_DIR = ana_dir
        
    def __str__(self):
        return "USER: \t {}\nICDATA_DIR: \t {}\nDATA_DIR: \t {}\nANA_DIR: \t {}\n".format(self.USER, self.ICDATA_DIR, self.DATA_DIR, self.ANA_DIR)