import os
import pickle


class Utils():
    def __init__(self,path):
        self.path=path

    def save_obj(self,obj, name):
        obj_path = os.path.join(self.path, name)
        # str=pickle.dumps(obj)
        with open(obj_path + '.pkl', 'wb') as f:
            pickle.dump(obj, f, pickle.DEFAULT_PROTOCOL)

    def load_obj(self,name):
        obj_path = os.path.join(self.path, name)
        with open(obj_path + '.pkl', 'rb') as f:
            return pickle.load(f)

    def write_to_csv(self, df, name):
        df.to_csv(name, index=False)

    def read_from_csv(self,name):
        csv_path=os.path.join(self.path, name)
        with open(csv_path, encoding="utf8") as csv_file:
            contents = csv_file.read()
        return contents
