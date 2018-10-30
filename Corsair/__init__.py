#!/usr/bin/env python3
"""init for the Worldbuilder module"""

import Worldbuilder.Load_objects as lo
from Worldbuilder.get_spell_levels import get_spell_levels
from Worldbuilder.Spell_cache.random_spell import random_spell


## classes
from Worldbuilder.Item_cache.Item import Item
from Worldbuilder.Item_cache.Armor import Armor
from Worldbuilder.Item_cache.Enchantment import Enchantment
from Worldbuilder.Item_cache.Potion import Potion
from Worldbuilder.Item_cache.Ring import Ring
from Worldbuilder.Item_cache.Rod import Rod
from Worldbuilder.Item_cache.Scroll import Scroll
from Worldbuilder.Item_cache.Staff import Staff
from Worldbuilder.Item_cache.Wand import Wand
from Worldbuilder.Item_cache.Weapon import Weapon
from Worldbuilder.Item_cache.Wonderous import Wonderous

from Worldbuilder.Location_cache.Location import Location
from Worldbuilder.Location_cache.Region import Region
from Worldbuilder.Location_cache.Shop import Shop
from Worldbuilder.Location_cache.Shopkeep import Shopkeep


## tables
from Worldbuilder.Roll_tables.random_magic_item import random_magic_item
from Worldbuilder.Roll_tables.armor_table import armor_table
from Worldbuilder.Roll_tables.potion_table import potion_table
from Worldbuilder.Roll_tables.ring_table import ring_table
from Worldbuilder.Roll_tables.rod_table import rod_table
from Worldbuilder.Roll_tables.scroll_table import scroll_table
from Worldbuilder.Roll_tables.staff_table import staff_table
from Worldbuilder.Roll_tables.wand_table import wand_table
from Worldbuilder.Roll_tables.weapon_table import weapon_table
from Worldbuilder.Roll_tables.wonderous_item_table import wonderous_table
from Worldbuilder.Roll_tables.store_list_table import store_list_table
from Worldbuilder.Roll_tables.items_for_shop_table import items_for_shop_table
from Worldbuilder.Roll_tables.item_spawn_chance_table import item_spawn_chance_table
from Worldbuilder.Roll_tables.store_by_loc_size import store_by_loc_size


## need an __all__ statement
# __all__ = ["echo", "surround", "reverse"]


## we need to initate some globals here
# SPELLS - a global dicitonary of spells to use
# SPELLS_BY_LEVEL - global dictionary of all the spells by level
## load spells
# try:
#     with open('Worldbuilder/pathfinder_files/spells.pkl', 'rb') as f:
#         SPELLS = pickle.load(f)
# except:
#     print('Could not load spells, loaded from table.\n')
#     SPELLS = lo.load_spells('Worldbuilder/pathfinder_files/spells.xls')
# SPELLS_BY_LEVEL = get_spell_levels(SPELLS)


# EVERYHTING -global dictionary that contains everything that has been made, stored by hashes
# SHOPS - all shops