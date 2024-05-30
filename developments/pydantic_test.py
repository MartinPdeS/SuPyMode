from pydantic import validate_call


@validate_call
def test(value: str):
    a = value
    print(a)


test(2.0)