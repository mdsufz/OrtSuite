from bs4 import BeautifulSoup as bs


def parser_ids(file):
    res = []
    with open(file, 'r') as f:
        lines = f.readlines()
    for line in lines:
        kID = line.strip().upper()
        if len(kID) > 0:
            res.append(kID)
    return res


def check_ec(ecs):
    ec_list = [x.lower() for x in ecs]
    if ec_list[0][0:2] == 'ec':
        ec_list = [x.replace('ec', '').strip() for x in ec_list]
    ec_list = [x.replace(' ', '') for x in ec_list]
    return ec_list


def make_blocks(url_list):
    # How many url each chunk should have
    n = 100
    res = [url_list[i * n:(i + 1) * n] for i in range((len(url_list) + n - 1) // n)]
    return res


# API links

def get_orthology_ids_url_from_map(pathway_id):
    URL = 'http://www.genome.jp'
    FUN = '/dbget-bin/get_linkdb?-t+orthology+path:'
    return URL + FUN + pathway_id


def get_gene_ids_url(orthology_id):
    URL = 'http://www.genome.jp'
    FUN = '/dbget-bin/get_linkdb?-t+genes+ko:'
    return URL + FUN + orthology_id


def get_orthology_url_from_ec(ec):
    URL = 'https://www.genome.jp'
    FUN = '/dbget-bin/get_linkdb?-t+orthology+ec:'  # ec:1.3.1.25'
    return URL + FUN + ec


def get_orthology_url_from_rn(rn):
    URL = 'https://www.genome.jp'
    FUN = '/dbget-bin/get_linkdb?-t+orthology+rn:'
    return URL + FUN + rn


def get_ko_url(ko):
    URL = 'https://www.genome.jp'
    FUN = '/dbget-bin/www_bget?ko:'
    return URL + FUN + ko


def get_ec_url(ec):
    URL = 'https://www.genome.jp'
    FUN = '/dbget-bin/www_bget?ec:'
    return URL + FUN + ec


def get_fastaProt_url(prot_id):
    URL = 'http://www.genome.jp'
    FUN = '/dbget-bin/www_bget?-f+-n+a+'
    return URL + FUN + prot_id


def get_fasta_url(gene_id):
    URL = 'http://www.genome.jp'
    FUN = '/dbget-bin/www_bget?-f+-n+n+'
    return URL + FUN + gene_id


def get_ec_url_from_ko(ko):
    URL = 'https://www.genome.jp'
    FUN = '/dbget-bin/get_linkdb?-t+enzyme+ko:'
    return URL + FUN + ko


def get_rn_url_from_ko(ko):
    URL = 'https://www.genome.jp'
    FUN = '/dbget-bin/get_linkdb?-t+rn+ko:'
    return URL + FUN + ko


def get_ec_url_from_rn(rn):
    URL = 'https://www.genome.jp'
    FUN = '/dbget-bin/get_linkdb?-t+enzyme+rn:'
    return URL + FUN + rn


def get_ids(response):
    try:
        html = response.text
        b = bs(html, features="html.parser")
        links = b.find_all('a')
        valid_link = lambda x: 'www_bget' in x.get('href')
        links = filter(valid_link, links)
        lista = [link.text for link in links]
        return lista
    except AttributeError:
        html = response.read()
        b = bs(html, features="html.parser")
        links = b.find_all('a')
        valid_link = lambda x: 'www_bget' in x.get('href')
        links = filter(valid_link, links)
        lista = [link.text for link in links]
        return lista
    else:
        return None


def get_fastaProt(response):
    try:
        html = bs(response.text, features="html.parser")
        return html.pre.text
    except:
        return None


def get_ko_name(response):
    html = response.read()
    b = bs(html, features="html.parser")
    rows = b.findAll("tr")
    for r in rows:
        lines = r.find("nobr")
        if lines:
            n = lines.text
            if n == 'Definition':
                definition = r.find("td").text
    return definition


def get_ec_names(response):
    names = []
    html = response.read()
    b = bs(html, features="html.parser")
    rows = b.findAll("tr")
    for r in rows:
        lines = r.find("nobr")
        if lines:
            n = lines.text
            if n == 'Name':
                cells = r.findAll("td")
                for cell in cells:
                    names = [x.strip() for x in cell.text.split(';')]
    return names
