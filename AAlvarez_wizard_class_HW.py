##RENAME: YourLastName_wizrd_class_HW.py
##Part of Python Project: 5 points
## 
## Part 1: Character class
##  Modify Attack so that it generates a (wil need to import random) random integer from 1 to self.attack

#self is used to to refer to the current object being called on

##  and subtracts that from other.hp

## Part 2: Wizard class
## Create a child class, Wizard that inherits the parent class Character
## Use super() to bring in all the data/methods from Character
## (1) Add a parameter magic to the constructor
## (2) Add a max_hp value in the constructor that is set to the original value
##     of the wizard hp
## Add 2 methods to wizard
## (1) heal - that makes hp = max_hp
## (2) fireball - It is like the Attack method, but subtracts a random
##     integer between 1 and self.magic to other.
import random

class Character:
    def __init__(self, name='', hp=0, attack=0):
        self.name = name
        self.hp = hp
        self.attack = attack

#Change operator overloader to print the attack value.    
    def __str__(self):
        return "Name: " + self.name + "\nHP:" + str(self.hp) + "\nAttack: " + str(self.attack)
        
#random.randint choose intergeer between 1 to attack value and damages other.hp  
    def Attack(self, other):
        damage = random.randint(1, self.attack) 
        other.hp-=damage 
        print(self.name +" attacked " +other.name +str(damage) + " damage points ")
        print(other.name +" has " + str(other.hp)+ " HP reamining  \n")
        
        
#part 2: creat a child class called wizard from parent class (Charachter)    
class Wizard(Character):
    def __init__(self, name='', hp=0, attack=0, magic=0):
        super().__init__(name, hp, attack)
        self.magic = magic
        self.max_hp = hp 
		
    def heal(self):
        self.hp = self.max_hp
        print(self.name + " has been healed to max HP ")

    def fireball(self, other):
        damage = random.randint(1, self.magic)
        other.hp -= damage
        print(self.name + " fireball attack on " + other.name + " for " + str(damage) + " damage points")
        print(other.name + " has " + str(other.hp) + " HP remaining. \n")
        
    def __str__(self):
        return super().__str__() + "\nMagic: " + str(self.magic) + "\nMax HP: " + str(self.max_hp)


# str will create a string of the characters qualities 
##TEST with the following code###
        

if __name__=="__main__":
    print("Ready...Set...Fight!")
    
    weezle=Wizard('Weezle',200,5,500)
    grognak=Character('Grognak',250,300)
    print(weezle, "\n")
    print(grognak)
    grognak.Attack(weezle)
    print(weezle, "\n")
    weezle.fireball(grognak)
    print(grognak, "\n")
    weezle.heal()
    weezle.fireball(grognak)
    print(weezle, "\n")
    print(grognak)


## OPTIONAL: For fun only - add a while loop to main that has them fight      
##           each other until one of their hp is less than 0.
##           It then declares the other the winner.
## print("Grognak crushes Weezle!") or print("Weezle vaporizes Grognak!")

    ##Set it so that it is a more equal fight
    weezle=Wizard('Weezle',2000,5,500)
    grognak=Character('Grognak',2500,450)


