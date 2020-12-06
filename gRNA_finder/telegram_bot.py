import logging

import telegram
from telegram.ext import MessageHandler, Filters, Updater, CommandHandler

from gRNA_finder.handlers import start, echo
from gRNA_finder.local import TOKEN

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

logger = logging.getLogger()
logger.setLevel(logging.INFO)

bot = telegram.Bot(token=TOKEN)
# print(bot.get_me())

updater = Updater(token=TOKEN, use_context=True)
dispatcher = updater.dispatcher

start_handler = CommandHandler('start', start)
dispatcher.add_handler(start_handler)


echo_handler = MessageHandler(Filters.text & (~Filters.command), echo)
dispatcher.add_handler(echo_handler)

updater.start_polling()
updater.stop()
# updater.idle()

# if __name__ == '__main__':
#     main()
