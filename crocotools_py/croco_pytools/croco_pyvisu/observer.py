class Subscriber:
    def __init__(self, who):
    	self.who = who
    def update(self, level, longitude, latitude):
        self.who.update(level, longitude, latitude)

        
class Publisher:
    def __init__(self):
        self.subscribers = dict()
    def register(self, who, callback=None):
        if callback == None:
            callback = getattr(who, 'update')
        self.subscribers[who] = callback
    def unregister(self, who):
        del self.subscribers[who]
    def dispatch(self, level, longitude, latitude):
        for subscriber, callback in self.subscribers.items():
            callback(level, longitude, latitude)
