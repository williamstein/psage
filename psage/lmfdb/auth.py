def userpass():
    import getpass
    username = raw_input('Mongodb username: ')
    password = getpass.getpass()
    return username, password
