import telebot
from telebot import types
import numpy as np
import os

from Bio import Entrez
import time



#handle = Entrez.efetch(db="nucleotide", id=7157, rettype="gb", retmode="text")

bot = telebot.TeleBot(os.environ.get('TOKEN'))

english_greetings = ['Hi', 'Hey', 'Hello', 'Good to see you', 'I greet you']
russian_greetings = ['Привет', 'Здравствуй', 'Приветствую', 'Салют', 'Рад тебя видеть']

#keyboard = telebot.types.ReplyKeyboardMarkup(False,True)
#keyboard.row('/start', '/sequence', '/help')'''
Entrez.email = 'nig.gov@gmail.com'
time.sleep(2)
@bot.message_handler(commands=['sequence'])
def sequence(message):
    ss = " ".join(message.text.split()[1:])
    if ss == "":
        bot.send_message(message.chat.id, "please enter a search term")
    else:
        bot.send_message(message.chat.id, "searching for '" + ss +"'")

        esh = Entrez.esearch(db = "gene", term = ss)
        esr = Entrez.read(esh)
        esh.close()
        #handle = Entrez.efetch(db="gene", id=7157, retmode="text")
        bot.send_message(message.chat.id, str(esr['Count']) + " results")
        if esr['Count'] == 0:
            bot.send_message(message.chat.id, "please try another term")
            return
        for i in range (3):
                top_id = str(esr['IdList'][i])
                handle = Entrez.efetch(db="sequences", id=str(top_id), rettype="gb", retmode="xml")
                data_dict = Entrez.read(handle)[0]
                handle.close()
                try:
                    bot.send_message(message.chat.id, "The result is " + str(data_dict['GBSeq_definition'] + " located in " + str(data_dict['GBSeq_locus']) + " of " + str(data_dict['GBSeq_organism']) + " genome"))
                    #bot.send_message(message.chat.id, "located in " + str(data_dict['GBSeq_locus']) + " of " + str(data_dict['GBSeq_organism']) + " genome")
                    bot.send_message(message.chat.id, "the " + str(data_dict['GBSeq_moltype']) + " sequence is :\n" + str(data_dict['GBSeq_sequence']))
                    #bot.send_message(message.chat.id, "comment: " + str(data_dict['GBSeq_comment']))
                except KeyError:
                    pass
time.sleep(2)
@bot.message_handler(commands=['geneid'])
def geneid(message):
    ss = " ".join(message.text.split()[1:])
    if ss == "":
        bot.send_message(message.chat.id, "please enter a search term")
    else:
        bot.send_message(message.chat.id, "searching for '" + ss +"'")
        esh = Entrez.esearch(db = "gene", term = ss)
        esr = Entrez.read(esh)
        esh.close()
        #handle = Entrez.efetch(db="gene", id=7157, retmode="text")
        bot.send_message(message.chat.id, str(esr['Count']) + " results")
        if esr['Count'] == 0:
            bot.send_message(message.chat.id, "please try another term")
            return
        for i in range (2):
                top_id = str(esr['IdList'][i])
                handle = Entrez.efetch(db="gene", id=str(top_id), rettype="gb", retmode="xml")
                data_dict = Entrez.read(handle)[0]
                handle.close()
                try:
                    bot.send_message(message.chat.id, 'The result is ' + str(data_dict["Entrezgene_gene"]['Gene-ref']['Gene-ref_desc']))
                    bot.send_message(message.chat.id, ' with gene ID ' +  str((data_dict['Entrezgene_track-info']['Gene-track']['Gene-track_geneid'])))
                    bot.send_message(message.chat.id, ' from organism ' + str(data_dict['Entrezgene_source']['BioSource']['BioSource_org']['Org-ref']['Org-ref_taxname']))
                except KeyError:
                    pass

@bot.message_handler(commands=['start'])
def send_welcome(message):
         bot.reply_to(message, f'''Привет, {message.from_user.first_name}. Инициализация прошла успешно. Для списка доступных команд напишите /help.
         Мои создатели ленивы и тяжелы на подъем, и я им под стать, так что если я чего-то не могу, то ''')


@bot.message_handler(commands=['help'])
def help_message(message):
    bot.send_message(message.chat.id, 'Как мной пользоваться: пишешь /sequence X , где X представляет некий ID последовательности на NCBI; пишешь /geneid X, где X - это ID гена.\nВ ответ я выдам тебе кое-какую информацию из первых найденных результатов.  ') #reply_markup=keyboard)


neponal = '''Excuse me?. Enter /sequence X or /geneid X, where X is your search item'''
@bot.message_handler(content_types=['text'])
def send_text(message):
        bot.send_message(message.chat.id, neponal)





'''    #bot.reply_to(message, message.text)
    for i in english_greetings:
        if message.text.capitalize() == i:
            bot.send_message(message.chat.id, english_greetings[int(np.random.randint(len(english_greetings)))], ',', {message.from_user.first_name})
        else:
            bot.send_message(message.chat.id, neponal)
            break
    for n in russian_greetings:
        if message.text.capitalize() == n:
            bot.send_message(message.chat.id, russian_greetings[int(np.random.randint(len(russian_greetings)))], ',', {message.from_user.first_name})''
        else:
            bot.send_message(message.chat.id, neponal)
            break'''



bot.polling(True)
