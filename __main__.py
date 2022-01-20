#!/usr/bin/env python3
import os
import time

from telebot import types
import telegram
import telegram.bot
from telegram.ext import Updater, Dispatcher, CallbackContext, CommandHandler, MessageHandler, Filters
from Bio import Entrez, SeqIO

import logging
import tempfile

# omfg = '''
# NAME ON NCBI      DATABASE FOR SEARCH
# BioProject        bioproject
# BioSample         biosample
# Biosystems        biosystems
# Books             books
# Conserved Domains cdd
# dbGaP             gap
# dbVar             dbvar
# Epigenomics       epigenomics
# EST               nucest
# Gene              gene
# Genome            genome
# GEO Datasets      gds
# GEO Profiles      geoprofiles
# GSS               nucgss
# HomoloGene        homologene
# MeSH              mesh
# NCBI C++ Toolkit  toolkit
# NCBI Web Site     ncbisearch
# NLM Catalog       nlmcatalog
# Nucleotide        nuccore
# OMIA              omia
# PopSet            popset
# Probe             probe
# Protein           protein
# Protein Clusters  proteinclusters
# PubChem BioAssay  pcassay
# PubChem Compound  pccompound
# PubChem Substance pcsubstance
# PubMed            pubmed
# PubMed Central    pmc
# SNP               snp
# SRA               sra
# Structure         structure
# Taxonomy          taxonomy
# UniGene           unigene
# UniSTS            unists
# '''



omg = '''
NAME ON NCBI      DATABASE FOR SEARCH
EST                                       nucest
GSS                                       nucgss
Nucleotide                            nuccore
Nucleotide                            nucleotide
PopSet (USE CAREFULLY!)  popset
Protein                                   protein
'''



ncbi_dbs = '''
nucest
nucgss
nuccore
popset
protein
nucleotide'''

def on_start(update: telegram.Update, context: CallbackContext) -> None:
    """ Greetz user, logs id
    """
    logging.log(level=logging.INFO,
                msg=f"@{update.effective_chat.username} started bot")
    bot: telegram.bot.Bot = context.bot
    bot.send_message(chat_id=update.effective_chat.id,
                     text="Welcome, to AVI bot! Made by the most AVId men!\nUsage: write /dbname searchitem, where dbname is the name of your database, and searchitem is your search item.\nReturns head of specified search query in genbank format from the selected database.")

commands_available = ["start", "help", "databases", "nucest", "nucgss", "nuccore", "popset", "protein", "nucleotide"]

def on_miss(update: telegram.Update, context: CallbackContext) -> None:
    bot = update.effective_chat.bot
    i = ''.join(update.message.text.split(' ')[0][1:])
    if i in commands_available:
            pass
    else:
            bot.send_message(chat_id=update.effective_chat.id,
                         text =   "Something's wrong. I do not get you. Enter /help if you need any")
            return

def on_db_available(update: telegram.Update, context: CallbackContext) -> None:
    bot = update.effective_chat.bot
    bot.send_message(chat_id=update.effective_chat.id,
                 text =  omg)

def on_help(update: telegram.Update, context: CallbackContext) -> None:
    bot = update.effective_chat.bot
    bot.send_message(chat_id=update.effective_chat.id,
                 text = 'To use me, just write /dbname searchitem, where dbname is the name of your database, and searchitem is your... search item.\nE.g. /protein p53 homo sapiens tumor supressor; /nucleotide betta splendens.\nTo get the list of databases available, enter /databases')


def on_any(update: telegram.Update, context: CallbackContext) -> None:
    """ Delivers sequences by Entrez query
    """
    sequences = update.message.text.split(' ')[1:]
    bot = update.effective_chat.bot
    logging.log(level=logging.INFO,
                msg=f"@{update.effective_chat.username} requested {''.join(update.message.text.split(' ')[0][1:])} [{' '.join(sequences)}]")
    # Fool proof
    if not len(sequences):
        bot.send_message(chat_id=update.effective_chat.id,
                         text="""Empty search, try harder!""")
        return

    if ''.join(update.message.text.split(' ')[0][1:]) in ncbi_dbs:
            pass
    else:
        bot.send_message(chat_id=update.effective_chat.id,
                         text="""Invalid database. Print /databases to get the list of databases available""")
        return
    # Making request
    handle = Entrez.read(Entrez.esearch(db=''.join(update.message.text.split(' ')[0][1:]),
                                        term=' '.join(sequences),
                                        retmode="xml"))
    bot.send_message(chat_id=update.effective_chat.id,
                     text=f"Found {handle['Count']} records, sending {len(handle['IdList'])}, " +
                          f"query as {handle['QueryTranslation']}")
    # Fetching
    logging.log(level=logging.INFO,
                msg=f"Fetching {len(handle['IdList'])} entities for @{update.effective_chat.username}")
    records = SeqIO.parse(Entrez.efetch(db=''.join(update.message.text.split(' ')[0][1:]),
                                        id=handle["IdList"],
                                        rettype="gb",
                                        retmode="text"),
                          "gb")
    # Sending
    files = []  # For further use
    with tempfile.TemporaryDirectory() as t_dir:
        for seq_record in records:
            path = os.path.join(t_dir, seq_record.id + ".gb")
            files.append(path)
            SeqIO.write(seq_record, path, "gb")
            with open(path, "rb") as payload:
                bot.send_document(chat_id=update.effective_chat.id,
                                  caption=f"{seq_record.id}\n{seq_record.name}\n{seq_record.description}",
                                  document=payload,
                                  disable_notification=True,
                                  reply_to_message_id=update.message.message_id)
            logging.log(level=logging.INFO,
                        msg=f"Sent {files[-1]} " + f"to @{update.effective_chat.username}")



if __name__ == "__main__":
    from os import environ
    logging.basicConfig(level=logging.INFO)

    updater = Updater(environ['TG_API_TOKEN'])

    Entrez.email = environ['ENTREZ_EMAIL']

    # handlers go here
    updater.dispatcher.add_handler(CommandHandler("start", on_start))
    updater.dispatcher.add_handler(CommandHandler("databases", on_db_available))
    updater.dispatcher.add_handler(CommandHandler("help", on_help))
    updater.dispatcher.add_handler(CommandHandler("protein", on_any))
    updater.dispatcher.add_handler(CommandHandler("nucleotide", on_any))
    updater.dispatcher.add_handler(CommandHandler("nuccore", on_any))
    updater.dispatcher.add_handler(CommandHandler("nucgss", on_any))
    updater.dispatcher.add_handler(CommandHandler("nucest", on_any))
    updater.dispatcher.add_handler(CommandHandler("popset", on_any))
    # updater.dispatcher.add_handler(CommandHandler("pubmed", on_any))
    # updater.dispatcher.add_handler(CommandHandler("structure", on_any))
    # updater.dispatcher.add_handler(CommandHandler("genome", on_any))
    # updater.dispatcher.add_handler(CommandHandler("books", on_any))
    # updater.dispatcher.add_handler(CommandHandler("cancerchromosomes", on_any))
    # updater.dispatcher.add_handler(CommandHandler("cdd", on_any))
    # updater.dispatcher.add_handler(CommandHandler("gap", on_any))
    # updater.dispatcher.add_handler(CommandHandler("domains", on_any))
    # updater.dispatcher.add_handler(CommandHandler("gene", on_any))
    # updater.dispatcher.add_handler(CommandHandler("genomeprj", on_any))
    # updater.dispatcher.add_handler(CommandHandler("gensat", on_any))
    # updater.dispatcher.add_handler(CommandHandler("geo", on_any))
    # updater.dispatcher.add_handler(CommandHandler("gds", on_any))
    # updater.dispatcher.add_handler(CommandHandler("homologene", on_any))
    # updater.dispatcher.add_handler(CommandHandler("journals", on_any))
    # updater.dispatcher.add_handler(CommandHandler("mesh", on_any))
    # updater.dispatcher.add_handler(CommandHandler("ncbisearch", on_any))
    # updater.dispatcher.add_handler(CommandHandler("nlmcatalog", on_any))
    # updater.dispatcher.add_handler(CommandHandler("omia", on_any))
    # updater.dispatcher.add_handler(CommandHandler("omim", on_any))
    # updater.dispatcher.add_handler(CommandHandler("pmc", on_any))
    # updater.dispatcher.add_handler(CommandHandler("probe", on_any))
    # updater.dispatcher.add_handler(CommandHandler("proteinclusters", on_any))
    # updater.dispatcher.add_handler(CommandHandler("pcassay", on_any))
    # updater.dispatcher.add_handler(CommandHandler("pccompound", on_any))
    # updater.dispatcher.add_handler(CommandHandler("pcsubstance", on_any))
    # updater.dispatcher.add_handler(CommandHandler("snp", on_any))
    # updater.dispatcher.add_handler(CommandHandler("taxonomy", on_any))
    # updater.dispatcher.add_handler(CommandHandler("toolkit", on_any))
    # updater.dispatcher.add_handler(CommandHandler("unigene", on_any))
    # updater.dispatcher.add_handler(CommandHandler("unists", on_any))


    # updater.dispatcher.add_handler(CommandHandler("help", on_help))
    # default
    updater.dispatcher.add_handler(MessageHandler(Filters.update, on_miss))

    updater.start_polling(poll_interval=2.5)
    updater.idle()
