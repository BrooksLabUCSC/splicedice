default_config = {
    "significance_threshold":0.05,
    "delta_threshold":0.05,
    "beta_exclude_01s":False,
    "":"",
    "":"",
    "":"",
    "":"",
    "":"",
    "":"",

}

def get_config(filename=None,extra_args=""):
    config = default_config
    if filename:
        with open(filename) as txt:
            for line in txt:
                row = line.rstrip().split('=')
                config[row[0]] = row[1]

    for item in extra_args.split(','):
        print(item)
        if item:
            attribute,value = (x.strip() for x in item.split('='))
            try:
                print("try")
                if '.' in value:
                    value = float(value)
                else:
                    value = int(value)
            except ValueError:
                pass
            config[attribute] = value
    return config
