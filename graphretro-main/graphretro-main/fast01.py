from fastapi import FastAPI
import uvicorn
from typing import Union  #  Union 类型，用于支持多种数据类型的参数注解
from pydantic import BaseModel

app = FastAPI()

class Item(BaseModel):  # 整合的输入要求    请求体
    name: str
    price: float
    is_offer: Union[bool, None] = None

@app.get("/")
def read_root():
    return {"Hello": "World"}

@app.get("/items/{item_id}")
def read_item(item_id: int, q: Union[str, None] = None):   # item_id变成int类型，q变成str类型或者None
    return {"item_id": item_id, "q": q}

@app.put("/items/{item_id}")
def update_item(item_id: int, item: Item):
    return {"item_name": item.name, "item_id": item_id}


if __name__ == '__main__':  # 如果本文件是一个启动文件
    uvicorn.run("fast01:app", port=8000,reload=True)  # 端口，debug开发模式，重写加载   http://127.0.0.1:8000/docs

