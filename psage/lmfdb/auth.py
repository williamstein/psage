from builtins import input
def userpass():
    import getpass
    username = input('Mongodb username: ')
    password = getpass.getpass()
    return username, password
