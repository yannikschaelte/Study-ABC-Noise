class Model(ABC):

    def __init__(self):
        pass

    def get_prior(self):
        pass

    def call(self):
        pass

    def call_noisy(self):
        pass

    def save_to_file(self):
        pass

    @staticmethod
    def load_from_file():
        pass
